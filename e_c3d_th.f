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
      subroutine e_c3d_th(co,nk,kon,lakonl,s,sm,
     &  ff,nelem,nmethod,rhcon,nrhcon,ielmat,ielorien,norien,orab,
     &  ntmat_,t0,t1,ithermal,vold,iperturb,nelemload,
     &  sideload,xload,nload,idist,iexpl,dtime,
     &  matname,mi,mass,stiffness,buckling,rhsi,intscheme,
     &  physcon,shcon,nshcon,cocon,ncocon,ttime,time,istep,iinc,
     &  xstiff,xloadold,reltime,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,springarea,
     &  plkcon,nplkcon,npmat_,ncmat_,elcon,nelcon,lakon,
     &  pslavsurf,pmastsurf,mortar,clearini,plicon,nplicon,ipkon,
     &  ielprop,prop,iponoel,inoel,sti,xstateini,xstate,nstate_,
     &  network,ipobody,xbody,ibody)
!
!     computation of the element matrix and rhs for the element with
!     the topology in konl
!
!     ff: rhs without temperature and eigenstress contribution
!
!     nmethod=0: check for positive Jacobian
!     nmethod=1: stiffness matrix + right hand side
!     nmethod=2: stiffness matrix + mass matrix
!     nmethod=3: static stiffness + buckling stiffness
!     nmethod=4: stiffness matrix + mass matrix
!
      implicit none
!
      integer mass,stiffness,buckling,rhsi
!
      character*8 lakonl,lakon(*)
      character*20 sideload(*)
      character*80 matname(*),amat
!
      integer konl(26),ifaceq(8,6),nelemload(2,*),nk,nelem,nmethod,
     &  mattyp,ithermal(*),iperturb(*),nload,idist,i,j,k,i1,j1,l,m,
     &  ii,jj,id,ipointer,ig,kk,istiff,iperm(20),ipompc(*),mi(*),
     &  nrhcon(*),ielmat(mi(3),*),ielorien(mi(3),*),nodempc(3,*),nmpc,
     &  ntmat_,nope,nopes,norien,iexpl,imat,mint2d,ikmpc(*),
     &  mint3d,ifacet(6,4),nopev,iorien,ilmpc(*),kode,jfaces,null,
     &  ifacew(8,5),intscheme,ipointeri,ipointerj,ncocon(2,*),
     &  nshcon(*),iinc,istep,jltyp,nfield,node,iflag,iscale,ielprop(*),
     &  nplkcon(0:ntmat_,*),nelcon(2,*),npmat_,ncmat_,i2,ipkon(*),
     &  iemchange,kon(*),mortar,nplicon(0:ntmat_,*),indexe,igauss,
     &  iponoel(*),inoel(2,*),nstate_,network,ipobody(2,*),ibody(3,*)
!
      real*8 co(3,*),xl(3,26),shp(4,26),xstiff(27,mi(1),*),
     &  s(60,60),w(3,3),ff(60),shpj(4,26),sinktemp,xs2(3,7),
     &  rhcon(0:1,ntmat_,*),dxsj2,temp,press,xloadold(2,*),
     &  orab(7,*),t0(*),t1(*),coords(3),c1,c2,reltime,prop(*),
     &  xl2(3,9),xsj2(3),shp2(7,9),vold(0:mi(2),*),xload(2,*),
     &  xi,et,ze,xsj,xsjj,sm(60,60),t1l,rho,summass,summ,ttime,time,
     &  sume,factorm,factore,alp,weight,pgauss(3),timeend(2),
     &  cocon(0:6,ntmat_,*),shcon(0:3,ntmat_,*),sph,coconloc(6),
     &  field,areaj,sax(60,60),ffax(60),coefmpc(*),tl2(8),
     &  voldl(0:mi(2),26),springarea(2,*),plkcon(0:2*npmat_,ntmat_,*),
     &  elcon(0:ncmat_,ntmat_,*),elconloc(ncmat_),pslavsurf(3,*),
     &  pmastsurf(2,*),clearini(3,9,*),plicon(0:2*npmat_,ntmat_,*),
     &  sti(6,mi(1),*),xstate(nstate_,mi(1),*),xbody(7,*),
     &  xstateini(nstate_,mi(1),*),heatnod,heatfac
!
      real*8 dtime,physcon(*)
!
!
!
      include "gauss.f"
!
      ifaceq=reshape((/4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/),(/8,6/))
      ifacet=reshape((/1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/),(/6,4/))
      ifacew=reshape((/1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/),(/8,5/))
      null=0
      iperm=(/5,6,7,8,1,2,3,4,13,14,15,16,9,10,11,12,17,18,19,20/)
!
      iflag=3
!
      timeend(1)=time
      timeend(2)=ttime+time
!
      summass=0.d0
      heatfac=0.d0
!
      indexe=ipkon(nelem)
      imat=ielmat(1,nelem)
      amat=matname(imat)
      if(norien.gt.0) then
         iorien=max(0,ielorien(1,nelem))
      else
         iorien=0
      endif
!
      if(lakonl(4:5).eq.'20') then
         nope=20
         nopev=8
         nopes=8
      elseif(lakonl(4:4).eq.'8') then
         nope=8
         nopev=8
         nopes=4
      elseif(lakonl(4:5).eq.'10') then
         nope=10
         nopev=4
         nopes=6
      elseif(lakonl(4:4).eq.'4') then
         nope=4
         nopev=4
         nopes=3
      elseif(lakonl(4:5).eq.'15') then
         nope=15
         nopev=6
      elseif(lakonl(4:4).eq.'6') then
         nope=6
         nopev=6
      elseif(lakonl(1:2).eq.'ES') then
         if(lakonl(7:7).eq.'C') then
            if(mortar.eq.0) then
               nope=ichar(lakonl(8:8))-47
               konl(nope+1)=kon(indexe+nope+1)
            elseif(mortar.eq.1) then
               nope=kon(indexe)
            endif
         else
            nope=ichar(lakonl(8:8))-47
         endif
      elseif((lakonl(1:2).eq.'D ').or.
     &       ((lakonl(1:1).eq.'D').and.(network.eq.1))) then
         nope=3
      endif
!
c      if(intscheme.eq.0) then
!
!        # of 2D and 3D integration points
!
         if(lakonl(4:5).eq.'8R') then
            mint2d=1
            mint3d=1
         elseif(lakonl(4:7).eq.'20RB') then
            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
               mint2d=4
               mint3d=50
            else
               mint2d=4
               call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &              null,xi,et,ze,weight)
            endif
         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
            if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'E')) then
               mint2d=2
               mint3d=4
            else
               mint2d=4
               mint3d=8
            endif
         elseif(lakonl(4:4).eq.'2') then
            mint2d=9
            mint3d=27
         elseif(lakonl(4:5).eq.'10') then
            mint2d=3
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            mint2d=1
            mint3d=1
         elseif(lakonl(4:5).eq.'15') then
            mint3d=9
         elseif(lakonl(4:4).eq.'6') then
            mint3d=2
         else
            mint3d=0
         endif
c      else
c!
c!        # of 3D integration points
c!
c         if((lakonl(4:4).eq.'8').or.(lakonl(4:4).eq.'2')) then
c            mint3d=27
c         elseif((lakonl(4:5).eq.'10').or.(lakonl(4:4).eq.'4')) then
c            mint3d=15
c         elseif((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')) then
c            mint3d=9
c         else
c            mint3d=0
c         endif
c!
c!        # of 2D integration points
c!
c         if(lakonl(4:5).eq.'8R') then
c            mint2d=1
c         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
c            if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'E')) then
c               mint2d=2
c            else
c               mint2d=4
c            endif
c         elseif(lakonl(4:4).eq.'2') then
c            mint2d=9
c         elseif(lakonl(4:5).eq.'10') then
c            mint2d=3
c         elseif(lakonl(4:4).eq.'4') then
c            mint2d=1
c         endif
c      endif
!
!     computation of the coordinates of the local nodes
!
      do i=1,nope
        konl(i)=kon(indexe+i)
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
!
!       initialisation for distributed forces
!
      if(rhsi.eq.1) then
        if(idist.ne.0) then
          do i=1,nope
            ff(i)=0.d0
          enddo
        endif
      endif
!
!     initialisation of sm
!
      if((mass.eq.1).or.(buckling.eq.1)) then
        do i=1,nope
          do j=i,nope
            sm(i,j)=0.d0
          enddo
        enddo
      endif
!
!     initialisation of s
!
      do i=1,nope
        do j=i,nope
          s(i,j)=0.d0
        enddo
      enddo
!
!     calculating the stiffness matrix for the contact spring elements
!
      if(mint3d.eq.0) then
         do i1=1,nope
            do i2=0,3
               voldl(i2,i1)=vold(i2,konl(i1))
            enddo
         enddo
!
         if(lakonl(1:2).eq.'ES') then
            if(lakonl(7:7).eq.'C') then
!
!              contact element
!
               kode=nelcon(1,imat)
               if(kode.eq.-51) then
                  if(mortar.eq.0) then
                     call springstiff_n2f_th(xl,voldl,s,imat,elcon,
     &                 nelcon,ncmat_,ntmat_,nope,kode,plkcon,nplkcon,
     &                 npmat_,iperturb,springarea(1,konl(nope+1)),mi,
     &                 timeend,matname,konl(nope),nelem,istep,iinc)
                  elseif(mortar.eq.1) then
                     jfaces=kon(indexe+nope+2)
                     igauss=kon(indexe+nope+1) 
                     node=0
                     call springstiff_f2f_th(xl,voldl,s,imat,elcon,
     &                    nelcon,ncmat_,ntmat_,nope,lakonl,kode,
     &                    elconloc,plicon,nplicon,npmat_,
     &                    springarea(1,igauss),
     &                    nmethod,mi,reltime,jfaces,igauss,pslavsurf,
     &                    pmastsurf,clearini,matname,plkcon,nplkcon,
     &                    node,nelem,istep,iinc,timeend)
                  endif
               endif
            elseif(lakonl(7:7).eq.'F') then
!
!              advective element
!
               call advecstiff(nope,voldl,ithermal,xl,nelemload,
     &              nelem,nload,lakon,xload,istep,time,ttime,dtime,
     &              sideload,vold,mi,xloadold,reltime,nmethod,s,
     &              iinc,iponoel,inoel,ielprop,prop,ielmat,shcon,
     &              nshcon,rhcon,nrhcon,ntmat_,ipkon,kon,cocon,ncocon,
     &              ipobody,xbody,ibody)
            endif
         elseif((lakonl(1:2).eq.'D ').or.
     &          ((lakonl(1:1).eq.'D').and.(network.eq.1))) then
            do i=1,3
               do j=1,i
                  s(i,j)=0.d0
               enddo
            enddo
            call networkstiff(voldl,s,imat,konl,mi,ntmat_,shcon,
     &           nshcon,rhcon,nrhcon)
         endif
         return
      endif
!
!     computation of the matrix: loop over the Gauss points
!
      do kk=1,mint3d
c         if(intscheme.eq.0) then
            if(lakonl(4:5).eq.'8R') then
               xi=gauss3d1(1,kk)
               et=gauss3d1(2,kk)
               ze=gauss3d1(3,kk)
               weight=weight3d1(kk)
            elseif(lakonl(4:7).eq.'20RB') then
               if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
                  xi=gauss3d13(1,kk)
                  et=gauss3d13(2,kk)
                  ze=gauss3d13(3,kk)
                  weight=weight3d13(kk)
               else
                  call beamintscheme(lakonl,mint3d,ielprop(nelem),prop,
     &                 kk,xi,et,ze,weight)
               endif
            elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &              then
               xi=gauss3d2(1,kk)
               et=gauss3d2(2,kk)
               ze=gauss3d2(3,kk)
               weight=weight3d2(kk)
            elseif(lakonl(4:4).eq.'2') then
               xi=gauss3d3(1,kk)
               et=gauss3d3(2,kk)
               ze=gauss3d3(3,kk)
               weight=weight3d3(kk)
            elseif(lakonl(4:5).eq.'10') then
               xi=gauss3d5(1,kk)
               et=gauss3d5(2,kk)
               ze=gauss3d5(3,kk)
               weight=weight3d5(kk)
            elseif(lakonl(4:4).eq.'4') then
               xi=gauss3d4(1,kk)
               et=gauss3d4(2,kk)
               ze=gauss3d4(3,kk)
               weight=weight3d4(kk)
            elseif(lakonl(4:5).eq.'15') then
               xi=gauss3d8(1,kk)
               et=gauss3d8(2,kk)
               ze=gauss3d8(3,kk)
               weight=weight3d8(kk)
            else
               xi=gauss3d7(1,kk)
               et=gauss3d7(2,kk)
               ze=gauss3d7(3,kk)
               weight=weight3d7(kk)
            endif
c         else
c            if((lakonl(4:4).eq.'8').or.(lakonl(4:4).eq.'2')) then
c               xi=gauss3d3(1,kk)
c               et=gauss3d3(2,kk)
c               ze=gauss3d3(3,kk)
c               weight=weight3d3(kk)
c            elseif((lakonl(4:5).eq.'10').or.(lakonl(4:4).eq.'4')) then
c               xi=gauss3d6(1,kk)
c               et=gauss3d6(2,kk)
c               ze=gauss3d6(3,kk)
c               weight=weight3d6(kk)
c            else
c               xi=gauss3d8(1,kk)
c               et=gauss3d8(2,kk)
c               ze=gauss3d8(3,kk)
c               weight=weight3d8(kk)
c            endif
c         endif
!
!           calculation of the shape functions and their derivatives
!           in the gauss point
!
         if(nope.eq.20) then
            if(lakonl(7:7).eq.'A') then
               call shape20h_ax(xi,et,ze,xl,xsj,shp,iflag)
            elseif((lakonl(7:7).eq.'E').or.(lakonl(7:7).eq.'S')) then
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
!           check the jacobian determinant
!
         if(xsj.lt.1.d-20) then
            write(*,*) '*ERROR in e_c3d_th: nonpositive jacobian'
            write(*,*) '       determinant in element',nelem
            write(*,*)
            xsj=dabs(xsj)
            nmethod=0
         endif
!
!           calculating the temperature in the integration
!           point
!
         t1l=0.d0
!
         if((lakonl(4:7).eq.'20RA').or.(lakonl(4:7).eq.'20RS').or.
     &      (lakonl(4:7).eq.'20RE')) then
            t1l=vold(0,konl(1))*(shp(4,1)+shp(4,5)+shp(4,17))
     &         +vold(0,konl(2))*(shp(4,2)+shp(4,6)+shp(4,18))
     &         +vold(0,konl(3))*(shp(4,3)+shp(4,7)+shp(4,19))
     &         +vold(0,konl(4))*(shp(4,4)+shp(4,8)+shp(4,20))
     &         +vold(0,konl(9))*(shp(4,9)+shp(4,13))
     &         +vold(0,konl(10))*(shp(4,10)+shp(4,14))
     &         +vold(0,konl(11))*(shp(4,11)+shp(4,15))
     &         +vold(0,konl(12))*(shp(4,12)+shp(4,16))
         else
            do i1=1,nope
               t1l=t1l+shp(4,i1)*vold(0,konl(i1))
            enddo
         endif
!
!           calculating the coordinates of the integration point
!           for material orientation purposes (for cylindrical
!           coordinate systems)
!
         if(iorien.gt.0) then
            do j=1,3
               pgauss(j)=0.d0
               do i1=1,nope
                  pgauss(j)=pgauss(j)+shp(4,i1)*xl(j,i1)
               enddo
            enddo
         endif
!
!           material data
!
         istiff=1
         call materialdata_th(cocon,ncocon,imat,iorien,pgauss,orab,
     &        ntmat_,coconloc,mattyp,t1l,rhcon,nrhcon,rho,shcon,
     &        nshcon,sph,xstiff,kk,nelem,istiff,mi(1))
!
!           incorporating the jacobian determinant in the shape
!           functions
!
         xsjj=dsqrt(xsj)
         do i1=1,nope
            shpj(1,i1)=shp(1,i1)*xsjj
            shpj(2,i1)=shp(2,i1)*xsjj
            shpj(3,i1)=shp(3,i1)*xsjj
            shpj(4,i1)=shp(4,i1)*xsj
         enddo
!
         c1=coconloc(1)*weight
         c2=rho*sph*weight
!
!           determination of the stiffness, and/or mass and/or
!           buckling matrix
!
         do jj=1,nope
!
            do ii=1,jj
!
!                   the following section calculates the static
!                   part of the stiffness matrix which, for buckling 
!                   calculations, is done in a preliminary static
!                   call
!
               if(mattyp.eq.1) then
!
                  s(ii,jj)=s(ii,jj)+c1*
     &                 (shpj(1,ii)*shpj(1,jj)+shpj(2,ii)*shpj(2,jj)
     &                 +shpj(3,ii)*shpj(3,jj))
!
               elseif(mattyp.eq.2) then
!
                  s(ii,jj)=s(ii,jj)+(coconloc(1)*shpj(1,ii)*shpj(1,jj)
     &                +coconloc(2)*shpj(2,ii)*shpj(2,jj)
     &                +coconloc(3)*shpj(3,ii)*shpj(3,jj))*weight
!
               else
!
                  do i1=1,3
                     do j1=1,3
                        w(i1,j1)=shpj(i1,ii)*shpj(j1,jj)
                     enddo
                  enddo
!
                  s(ii,jj)=s(ii,jj)+
     &                (coconloc(1)*w(1,1)+
     &                 coconloc(4)*(w(1,2)+w(2,1))+
     &                 coconloc(2)*w(2,2)+
     &                 coconloc(5)*(w(1,3)+w(3,1))+
     &                 coconloc(6)*(w(2,3)+w(3,2))+
     &                 coconloc(3)*w(3,3))*weight
!
               endif
!
!                     mass matrix
!
               if(mass.eq.1) then
                  sm(ii,jj)=sm(ii,jj)
     &                 +c2*shpj(4,ii)*shp(4,jj)
               endif
!
            enddo
         enddo
!
!           computation of the right hand side
!
         if(rhsi.eq.1) then
!
!           distributed heat flux
!
            if(nload.gt.0) then
               call nident2(nelemload,nelem,nload,id)
               areaj=xsj*weight
               do
                  if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
                  if(sideload(id)(1:2).ne.'BF') then
                     id=id-1
                     cycle
                  endif
                  if(sideload(id)(3:4).eq.'NU') then
                     if(ithermal(2).eq.2) then
                        do j=1,3
                           pgauss(j)=0.d0
                           do i1=1,nope
                              pgauss(j)=pgauss(j)+
     &                             shp(4,i1)*xl(j,i1)
                           enddo
                        enddo
                     else
                        do j=1,3
                           pgauss(j)=0.d0
                           do i1=1,nope
                              pgauss(j)=pgauss(j)+
     &                             shp(4,i1)*(xl(j,i1)+vold(j,konl(i1)))
                           enddo
                        enddo
                     endif
                     jltyp=1
                     iscale=1
                     call dflux(xload(1,id),t1l,istep,iinc,timeend,
     &                 nelem,kk,pgauss,jltyp,temp,press,sideload(id),
     &                 areaj,vold,co,lakonl,konl,ipompc,nodempc,coefmpc,
     &                 nmpc,ikmpc,ilmpc,iscale,mi,sti,xstateini,xstate,
     &                 nstate_,dtime)
                     if((nmethod.eq.1).and.(iscale.ne.0))
     &                    xload(1,id)=xloadold(1,id)+
     &                   (xload(1,id)-xloadold(1,id))*reltime
                  endif
                  do jj=1,nope
                     ff(jj)=ff(jj)+xload(1,id)*shpj(4,jj)*weight
                  enddo
                  exit
               enddo
            endif
         endif
!
      enddo
!
      if((buckling.eq.0).and.(nload.ne.0)) then
         iflag=2
!
!       distributed loads
!
         call nident2(nelemload,nelem,nload,id)
         do
            if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
            if((sideload(id)(1:1).ne.'F').and.
     &         (sideload(id)(1:1).ne.'R').and.
     &         (sideload(id)(1:1).ne.'S')) then
               id=id-1
               cycle
            endif
c            read(sideload(id)(2:2),'(i1)') ig
            ig=ichar(sideload(id)(2:2))-48
!
!
!         treatment of wedge faces
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
          if((nope.eq.20).or.(nope.eq.8)) then
             do i=1,nopes
                tl2(i)=vold(0,konl(ifaceq(i,ig)))
                if(ithermal(2).eq.2) then
                   do j=1,3
                      xl2(j,i)=co(j,konl(ifaceq(i,ig)))
                   enddo
                else
                   do j=1,3
                      xl2(j,i)=co(j,konl(ifaceq(i,ig)))+
     &                     vold(j,konl(ifaceq(i,ig)))
                   enddo
                endif
             enddo
          elseif((nope.eq.10).or.(nope.eq.4)) then
             do i=1,nopes
                tl2(i)=vold(0,konl(ifacet(i,ig)))
                if(ithermal(2).eq.2) then
                   do j=1,3
                      xl2(j,i)=co(j,konl(ifacet(i,ig)))
                   enddo
                else
                   do j=1,3
                      xl2(j,i)=co(j,konl(ifacet(i,ig)))+
     &                     vold(j,konl(ifacet(i,ig)))
                   enddo
                endif
             enddo
          else
             do i=1,nopes
                tl2(i)=vold(0,konl(ifacew(i,ig)))
                if(ithermal(2).eq.2) then
                   do j=1,3
                      xl2(j,i)=co(j,konl(ifacew(i,ig)))
                   enddo
                else
                   do j=1,3
                      xl2(j,i)=co(j,konl(ifacew(i,ig)))+
     &                     vold(j,konl(ifacew(i,ig)))
                   enddo
                endif
             enddo
          endif
!
!         storing envnode values into xload for non-cavity
!         radiation
!
          if((sideload(id)(1:1).eq.'R').and.
     &       (sideload(id)(3:4).ne.'CR').and.
     &       (nelemload(2,id).gt.0)) then
             xload(2,id)=vold(0,nelemload(2,id))-physcon(1)
          endif
!
          do i=1,mint2d
!
!            copying the sink temperature to ensure the same
!            value in each integration point (sinktemp can be
!            changed in subroutine film: requirement from the
!            thermal people)
!
             sinktemp=xload(2,id)
!
             if((lakonl(4:5).eq.'8R').or.
     &            ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
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
     &               ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                xi=gauss2d5(1,i)
                et=gauss2d5(2,i)
                weight=weight2d5(i)
             elseif((lakonl(4:4).eq.'4').or.
     &               ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                xi=gauss2d4(1,i)
                et=gauss2d4(2,i)
                weight=weight2d4(i)
             endif
!
c             if(nopes.eq.9) then
c                call shape9q(xi,et,xl2,xsj2,xs2,shp2,iflag)
             if(nopes.eq.8) then
                call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
             elseif(nopes.eq.4) then
                call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
             elseif(nopes.eq.6) then
                call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
c             elseif(nopes.eq.7) then
c                call shape7tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
             else
                call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
             endif
!
             dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &            xsj2(3)*xsj2(3))
             areaj=dxsj2*weight
!
             temp=0.d0
             do j=1,nopes
                temp=temp+tl2(j)*shp2(4,j)
             enddo
!
!            for nonuniform load: determine the coordinates of the
!            point (transferred into the user subroutine)
!
             if((sideload(id)(3:4).eq.'NU').or.
     &          (sideload(id)(5:6).eq.'NU')) then
                do k=1,3
                   coords(k)=0.d0
                   do j=1,nopes
                      coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                   enddo
                enddo
                jltyp=ichar(sideload(id)(2:2))-48
                jltyp=jltyp+10
                if(sideload(id)(1:1).eq.'S') then
                   iscale=1
                   call dflux(xload(1,id),temp,istep,iinc,timeend,
     &               nelem,i,coords,jltyp,temp,press,sideload(id),
     &               areaj,vold,co,lakonl,konl,
     &               ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,mi,
     &               sti,xstateini,xstate,nstate_,dtime)
                elseif(sideload(id)(1:1).eq.'F') then
                   node=nelemload(2,id)
                   call film(xload(1,id),sinktemp,temp,istep,
     &               iinc,timeend,nelem,i,coords,jltyp,field,nfield,
     &               sideload(id),node,areaj,vold,mi,
     &               ipkon,kon,lakon,iponoel,inoel,ielprop,prop,ielmat,
     &               shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon,
     &               ipobody,xbody,ibody,heatnod,heatfac)
                elseif(sideload(id)(1:1).eq.'R') then
                   call radiate(xload(1,id),xload(2,id),temp,istep,
     &               iinc,timeend,nelem,i,coords,jltyp,field,nfield,
     &               sideload(id),node,areaj,vold,mi,iemchange)
                endif
             endif
                
!
             do k=1,nopes
                if((nope.eq.20).or.(nope.eq.8)) then
                   ipointer=ifaceq(k,ig)
                elseif((nope.eq.10).or.(nope.eq.4)) then
                   ipointer=ifacet(k,ig)
                else
                   ipointer=ifacew(k,ig)
                endif
                if(sideload(id)(1:1).eq.'S') then
!
!        flux INTO the face is positive (input deck convention)
!        this is different from the convention in the theory
!
                   ff(ipointer)=ff(ipointer)+shp2(4,k)*xload(1,id)
     &                  *areaj
                elseif(sideload(id)(1:1).eq.'F') then
                   ff(ipointer)=ff(ipointer)+shp2(4,k)*areaj*
     &                  (heatfac-xload(1,id)*(temp-sinktemp))
                elseif(sideload(id)(1:1).eq.'R') then
                   ff(ipointer)=ff(ipointer)-shp2(4,k)*physcon(2)*
     &                  xload(1,id)*((temp-physcon(1))**4-
     &                  (xload(2,id)-physcon(1))**4)*
     &                  areaj
                endif
             enddo
!
             do ii=1,nopes
                if((nope.eq.20).or.(nope.eq.8)) then
                   ipointeri=ifaceq(ii,ig)
                elseif((nope.eq.10).or.(nope.eq.4)) then
                   ipointeri=ifacet(ii,ig)
                else
                   ipointeri=ifacew(ii,ig)
                endif
                do jj=1,nopes
                   if((nope.eq.20).or.(nope.eq.8)) then
                      ipointerj=ifaceq(jj,ig)
                   elseif((nope.eq.10).or.(nope.eq.4)) then
                      ipointerj=ifacet(jj,ig)
                   else
                      ipointerj=ifacew(jj,ig)
                   endif
                   if(ipointeri.gt.ipointerj) cycle
                   if(sideload(id)(1:1).eq.'F') then
                      s(ipointeri,ipointerj)=s(ipointeri,ipointerj)+
     &                  shp2(4,ii)*shp2(4,jj)*xload(1,id)*areaj
                   elseif(sideload(id)(1:1).eq.'R') then
                      s(ipointeri,ipointerj)=s(ipointeri,ipointerj)+
     &                     shp2(4,ii)*shp2(4,jj)*physcon(2)*xload(1,id)*
     &                     areaj*4.d0*(temp-physcon(1))**3
                   endif
                enddo
             enddo
!
          enddo
!
          id=id-1
       enddo
      endif
!
!     for axially symmetric and plane stress/strain elements: 
!     complete s and sm
!
c      if(intscheme.eq.0) then
         if(((lakonl(4:5).eq.'8 ').or.
     &        ((lakonl(4:6).eq.'20R').and.(lakonl(7:8).ne.'BR'))).and.
     &        ((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'E'))) then
            do i=1,20
               do j=i,20
                  k=iperm(i)
                  l=iperm(j)
                  if(k.gt.l) then
                     m=k
                     k=l
                     l=m
                  endif
                  sax(i,j)=s(k,l)
               enddo
            enddo
            do i=1,20
               do j=i,20
                  s(i,j)=s(i,j)+sax(i,j)
               enddo
            enddo
!     
!     special treatment of plane stress elements since lateral
!     heating is allowed (orthogonal to the plane)
!     
            if(nload.ne.0) then
               do i=1,20
                  k=iperm(i)
                  ffax(i)=ff(k)
               enddo
               do i=1,20
                  ff(i)=ff(i)+ffax(i)
               enddo
            endif
!     
            if(mass.eq.1) then
               do i=1,20
                  do j=i,20
                     k=iperm(i)
                     l=iperm(j)
                     if(k.gt.l) then
                        m=k
                        k=l
                        l=m
                     endif
                     sax(i,j)=sm(k,l)
                  enddo
               enddo
               do i=1,20
                  do j=i,20
                     sm(i,j)=sm(i,j)+sax(i,j)
                  enddo
               enddo
            endif
         endif
c      else
c!
c!        only 2-D scheme is reduced for intscheme=1
c!
c         if(((lakonl(4:5).eq.'8 ').or.
c     &        ((lakonl(4:6).eq.'20R').and.(lakonl(7:8).ne.'BR'))).and.
c     &        ((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'E'))) then
c!     
c!           special treatment of plane stress elements since lateral
c!           heating is allowed (orthogonal to the plane)
c!     
c            if(nload.ne.0) then
c               do i=1,20
c                  k=iperm(i)
c                  ffax(i)=ff(k)
c               enddo
c               do i=1,20
c                  ff(i)=ff(i)+ffax(i)
c               enddo
c            endif
c         endif
c      endif
!
      if((mass.eq.1).and.(iexpl.gt.1)) then
!
!        scaling the diagonal terms of the mass matrix such that the total
!        mass is right (LUMPING; for explicit dynamic calculations)
!
         sume=0.d0
         summ=0.d0
         do i=1,nopev
            sume=sume+sm(i,i)
         enddo
         do i=nopev+1,nope
            summ=summ+sm(i,i)
         enddo
!
         if(nope.eq.20) then
            alp=.2917d0
         elseif(nope.eq.10) then
            alp=0.1203d0
         elseif(nope.eq.15) then
            alp=0.2141d0
         endif
!
         if((nope.eq.20).or.(nope.eq.10).or.
     &      (nope.eq.15)) then
            factore=summass*alp/(1.d0+alp)/sume
            factorm=summass/(1.d0+alp)/summ
         else
            factore=summass/sume
         endif
!
         do i=1,nopev
            sm(i,i)=sm(i,i)*factore
         enddo
         do i=nopev+1,nope
            sm(i,i)=sm(i,i)*factorm
         enddo
!
      endif
!
      return
      end



