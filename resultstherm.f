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
      subroutine resultstherm(co,kon,ipkon,lakon,v,
     &  elcon,nelcon,rhcon,nrhcon,ielmat,ielorien,norien,orab,
     &  ntmat_,t0,iperturb,fn,shcon,nshcon,
     &  iout,qa,vold,ipompc,nodempc,coefmpc,nmpc,
     &  dtime,time,ttime,plkcon,nplkcon,xstateini,xstiff,xstate,npmat_,
     &  matname,mi,ncmat_,nstate_,cocon,ncocon,
     &  qfx,ikmpc,ilmpc,istep,iinc,springarea,
     &  calcul_fn,calcul_qa,nal,nea,neb,ithermal,nelemload,nload,
     &  nmethod,reltime,sideload,xload,xloadold,pslavsurf,
     &  pmastsurf,mortar,clearini,plicon,nplicon,ielprop,prop,
     &  iponoel,inoel,network,ipobody,xbody,ibody)
!
!     calculates the heat flux and the material tangent at the integration
!     points and the internal concentrated flux at the nodes
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*20 sideload(*)
      character*80 amat,matname(*)
!
      integer kon(*),konl(26),iperm(20),ikmpc(*),ilmpc(*),mi(*),
     &  nelcon(2,*),nrhcon(*),ielmat(mi(3),*),ielorien(mi(3),*),
     &  ntmat_,ipkon(*),ipompc(*),nodempc(3,*),mortar,igauss,
     &  ncocon(2,*),iflag,nshcon(*),istep,iinc,mt,mattyp,
     &  i,j,k,m1,kk,i1,m3,indexe,nope,norien,iperturb(*),iout,
     &  nal,nmpc,kode,imat,mint3d,iorien,istiff,ncmat_,nstate_,
     &  nplkcon(0:ntmat_,*),npmat_,calcul_fn,calcul_qa,nea,neb,
     &  nelemload(2,*),nload,ithermal(*),nmethod,nopered,
     &  jfaces,node,nplicon(0:ntmat_,*),null,ielprop(*),
     &  iponoel(*),inoel(2,*),network,ipobody(2,*),ibody(3,*)
!
      real*8 co(3,*),v(0:mi(2),*),shp(4,26),reltime,
     &  xl(3,26),vl(0:mi(2),26),elcon(0:ncmat_,ntmat_,*),
     &  rhcon(0:1,ntmat_,*),qfx(3,mi(1),*),orab(7,*),
     &  rho,fn(0:mi(2),*),tnl(19),timeend(2),q(0:mi(2),26),
     &  vkl(0:3,3),t0(*),vold(0:mi(2),*),coefmpc(*),
     &  springarea(2,*),elconloc(ncmat_),cocon(0:6,ntmat_,*),
     &  shcon(0:3,ntmat_,*),sph,c1,xi,et,ze,xsj,qa(*),t0l,t1l,dtime,
     &  weight,pgauss(3),coconloc(6),qflux(3),time,ttime,
     &  t1lold,plkcon(0:2*npmat_,ntmat_,*),xstiff(27,mi(1),*),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),
     &  xload(2,*),xloadold(2,*),clearini(3,9,*),pslavsurf(3,*),
     &  pmastsurf(6,*),plicon(0:2*npmat_,ntmat_,*),prop(*),
     &  xbody(7,*)
!
      include "gauss.f"
!
      iflag=3
      null=0
      iperm=(/5,6,7,8,1,2,3,4,13,14,15,16,9,10,11,12,17,18,19,20/)
!
      mt=mi(2)+1
!
!     calculation of temperatures and thermal flux
!
      nal=0
!
      do i=nea,neb
!
!        only structural elements (no fluid elements)
!
         if((ipkon(i).lt.0).or.(lakon(i)(1:1).eq.'F')) cycle
!
         if(lakon(i)(1:7).eq.'DCOUP3D') cycle
!
!        strainless reactivated elements are labeled by a negative
!        value of ielmat. If the step is mechanical or thermo-mechanical
!        resultsmech is called before and the negative value is
!        reversed. If the step is purely thermal the negative sign
!        persists and has to be reverted here.
!
         ielmat(1,i)=abs(ielmat(1,i))
!
         imat=ielmat(1,i)
         amat=matname(imat)
         if(norien.gt.0) then
            iorien=max(0,ielorien(1,i))
         else
            iorien=0
         endif
!
         indexe=ipkon(i)
         if(lakon(i)(4:5).eq.'20') then
            nope=20
         elseif(lakon(i)(4:4).eq.'8') then
            nope=8
         elseif(lakon(i)(4:5).eq.'10') then
            nope=10
         elseif(lakon(i)(4:4).eq.'4') then
            nope=4
         elseif(lakon(i)(4:5).eq.'15') then
            nope=15
         elseif(lakon(i)(4:4).eq.'6') then
            nope=6
         elseif((lakon(i)(1:2).eq.'ES').and.(lakon(i)(7:7).ne.'A')) then
!
!           contact spring and advection elements (no dashpot elements
!           = ED... elements and no genuine spring elements)
!
            if(lakon(i)(7:7).eq.'C') then
!
!              contact spring elements
!
               if(mortar.eq.1) then
!
!                 face-to-face penalty
!
                  nope=kon(ipkon(i))
               elseif(mortar.eq.0) then
!
!                 node-to-face penalty
!
                  nope=ichar(lakon(i)(8:8))-47
                  konl(nope+1)=kon(indexe+nope+1)
               endif
            else
!
!              advection elements
!
               nope=ichar(lakon(i)(8:8))-47
            endif
c            nope=ichar(lakon(i)(8:8))-47
c!
c!           local contact spring number
c!
c            if(lakon(i)(7:7).eq.'C') konl(nope+1)=kon(indexe+nope+1)
         elseif((lakon(i)(1:2).eq.'D ').or.
     &          ((lakon(i)(1:1).eq.'D').and.(network.eq.1))) then
!
!           no entry or exit elements
!
            if((kon(indexe+1).eq.0).or.(kon(indexe+3).eq.0)) cycle
            nope=3
         else
            cycle
         endif
!
         if(lakon(i)(4:5).eq.'8R') then
            mint3d=1
         elseif(lakon(i)(4:7).eq.'20RB') then
            if((lakon(i)(8:8).eq.'R').or.(lakon(i)(8:8).eq.'C')) then
               mint3d=50
            else
               call beamintscheme(lakon(i),mint3d,ielprop(i),prop,
     &              null,xi,et,ze,weight)
            endif
         elseif((lakon(i)(4:4).eq.'8').or.
     &          (lakon(i)(4:6).eq.'20R')) then
            if(lakon(i)(6:7).eq.'RA') then
               mint3d=4
            else
               mint3d=8
            endif
         elseif(lakon(i)(4:4).eq.'2') then
            mint3d=27
         elseif(lakon(i)(4:5).eq.'10') then
            mint3d=4
         elseif(lakon(i)(4:4).eq.'4') then
            mint3d=1
         elseif(lakon(i)(4:5).eq.'15') then
            mint3d=9
         elseif(lakon(i)(4:4).eq.'6') then
            mint3d=2
         else
            mint3d=0
         endif
!
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               xl(k,j)=co(k,konl(j))
               vl(k,j)=v(k,konl(j))
            enddo
            vl(0,j)=v(0,konl(j))
         enddo
!
!        q contains the nodal forces per element; initialisation of q
!
         if((iperturb(1).ge.2).or.((iperturb(1).le.0).and.(iout.lt.1))) 
     &          then
            do m1=1,nope
               q(0,m1)=fn(0,konl(m1))
            enddo
         endif
!
!        calculating the concentrated flux for the contact elements
!
         if(mint3d.eq.0) then
!
            lakonl=lakon(i)
!
!           initialization of tnl
!
            do j=1,nope
               tnl(j)=0.d0
            enddo
!
!           spring elements (including contact springs)
!     
            if(lakonl(2:2).eq.'S') then
               if(lakonl(7:7).eq.'C') then
!
!                 contact element
!
                  kode=nelcon(1,imat)
                  if(kode.eq.-51) then
                     timeend(1)=time
                     timeend(2)=ttime+time
                     if(mortar.eq.0) then
                        call springforc_n2f_th(xl,vl,imat,elcon,nelcon,
     &                    tnl,ncmat_,ntmat_,nope,kode,elconloc,
     &                    plkcon,nplkcon,npmat_,mi,
     &                    springarea(1,konl(nope+1)),timeend,matname,
     &                    konl(nope),i,istep,iinc,iperturb)
                     elseif(mortar.eq.1) then
                        jfaces=kon(indexe+nope+2)
                        igauss=kon(indexe+nope+1)
                        node=0
                        call springforc_f2f_th(xl,vl,imat,elcon,nelcon,
     &                    tnl,ncmat_,ntmat_,nope,lakonl,kode,elconloc,
     &                    plicon,nplicon,npmat_,mi,springarea(1,igauss),
     &                    nmethod,reltime,jfaces,igauss,
     &                    pslavsurf,pmastsurf,clearini,timeend,istep,
     &                    iinc,plkcon,nplkcon,node,i,matname)
                     endif
                  endif
               elseif(lakonl(7:7).eq.'F') then
!
!                 advective element
!
                  call advecforc(nope,vl,ithermal,xl,nelemload,
     &                 i,nload,lakon,xload,istep,time,ttime,
     &                 dtime,sideload,v,mi,xloadold,reltime,nmethod,
     &                 tnl,iinc,iponoel,inoel,ielprop,prop,ielmat,shcon,
     &                 nshcon,rhcon,nrhcon,ntmat_,ipkon,kon,cocon,
     &                 ncocon,ipobody,xbody,ibody)
               endif
!
            elseif((lakonl(1:2).eq.'D ').or.
     &             ((lakonl(1:1).eq.'D').and.(network.eq.1))) then
!
!              generic networkelement
!
               call networkforc(vl,tnl,imat,konl,mi,ntmat_,shcon,
     &              nshcon,rhcon,nrhcon)
               
            endif
!
            do j=1,nope
               fn(0,konl(j))=fn(0,konl(j))+tnl(j)
            enddo
         endif
!
         do kk=1,mint3d
            if(lakon(i)(4:5).eq.'8R') then
               xi=gauss3d1(1,kk)
               et=gauss3d1(2,kk)
               ze=gauss3d1(3,kk)
               weight=weight3d1(kk)
            elseif(lakon(i)(4:7).eq.'20RB') then
               if((lakon(i)(8:8).eq.'R').or.(lakon(i)(8:8).eq.'C')) then
                  xi=gauss3d13(1,kk)
                  et=gauss3d13(2,kk)
                  ze=gauss3d13(3,kk)
                  weight=weight3d13(kk)
               else
                  call beamintscheme(lakon(i),mint3d,ielprop(i),
     &                 prop,kk,xi,et,ze,weight)
               endif
            elseif((lakon(i)(4:4).eq.'8').or.
     &             (lakon(i)(4:6).eq.'20R'))
     &        then
               xi=gauss3d2(1,kk)
               et=gauss3d2(2,kk)
               ze=gauss3d2(3,kk)
               weight=weight3d2(kk)
            elseif(lakon(i)(4:4).eq.'2') then
               xi=gauss3d3(1,kk)
               et=gauss3d3(2,kk)
               ze=gauss3d3(3,kk)
               weight=weight3d3(kk)
            elseif(lakon(i)(4:5).eq.'10') then
               xi=gauss3d5(1,kk)
               et=gauss3d5(2,kk)
               ze=gauss3d5(3,kk)
               weight=weight3d5(kk)
            elseif(lakon(i)(4:4).eq.'4') then
               xi=gauss3d4(1,kk)
               et=gauss3d4(2,kk)
               ze=gauss3d4(3,kk)
               weight=weight3d4(kk)
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
            c1=xsj*weight
!
!                 vkl(m2,m3) contains the derivative of the m2-
!                 component of the displacement with respect to
!                 direction m3
!
            do m3=1,3
               vkl(0,m3)=0.d0
            enddo
!
            do m1=1,nope
               do m3=1,3
                  vkl(0,m3)=vkl(0,m3)+shp(m3,m1)*vl(0,m1)
               enddo
            enddo
!
            kode=ncocon(1,imat)
!
!              calculating the temperature difference in
!              the integration point
!
            t1lold=0.d0
            t1l=0.d0
            if((lakon(i)(4:5).eq.'8 ').or.
     &         (lakon(i)(4:5).eq.'8I')) then
               do i1=1,8
                  t1lold=t1lold+vold(0,konl(i1))/8.d0
                  t1l=t1l+v(0,konl(i1))/8.d0
               enddo
            elseif(lakon(i)(4:6).eq.'20 ') then
               nopered=20
               call lintemp_th1(vold,konl,nopered,kk,t1lold,mi)
               call lintemp_th1(v,konl,nopered,kk,t1l,mi)
            elseif(lakon(i)(4:6).eq.'10T') then
               call linscal10(vold,konl,t1lold,mi(2),shp)
               call linscal10(v,konl,t1l,mi(2),shp)
            else
               do i1=1,nope
                  t1lold=t1lold+shp(4,i1)*vold(0,konl(i1))
                  t1l=t1l+shp(4,i1)*v(0,konl(i1))
               enddo
            endif
!
!           calculating the coordinates of the integration point
!           for material orientation purposes (for cylindrical
!           coordinate systems)
!
            if((iorien.gt.0).or.(kode.le.-100)) then
               do j=1,3
                  pgauss(j)=0.d0
                  do i1=1,nope
                     pgauss(j)=pgauss(j)+shp(4,i1)*co(j,konl(i1))
                  enddo
               enddo
            endif
!
!                 material data; for linear elastic materials
!                 this includes the calculation of the stiffness
!                 matrix
!
            istiff=0
!
            call materialdata_th(cocon,ncocon,imat,iorien,pgauss,orab,
     &           ntmat_,coconloc,mattyp,t1l,rhcon,nrhcon,rho,shcon,
     &           nshcon,sph,xstiff,kk,i,istiff,mi(1))
!
            call thermmodel(amat,i,kk,kode,coconloc,vkl,dtime,
     &           time,ttime,mi(1),nstate_,xstateini,xstate,qflux,xstiff,
     &           iorien,pgauss,orab,t1l,t1lold,vold,co,lakon(i),konl,
     &           ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,nmethod,
     &           iperturb)
! 
            qfx(1,kk,i)=qflux(1)
            qfx(2,kk,i)=qflux(2)
            qfx(3,kk,i)=qflux(3)
            if(lakon(i)(6:7).eq.'RA') then
               qfx(1,kk+4,i)=qflux(1)
               qfx(2,kk+4,i)=qflux(2)
               qfx(3,kk+4,i)=qflux(3)
            endif
!
!           calculation of the nodal flux
!
            if(calcul_fn.eq.1)then
!
!                    calculating fn using skl
!
               if(lakon(i)(6:7).eq.'RA') then
                  do m1=1,nope
                     fn(0,konl(m1))=fn(0,konl(m1))
     &                  -c1*(qflux(1)*(shp(1,m1)+shp(1,iperm(m1)))
     &                      +qflux(2)*(shp(2,m1)+shp(2,iperm(m1)))
     &                      +qflux(3)*(shp(3,m1)+shp(3,iperm(m1))))
                  enddo
               else
                  do m1=1,nope
                     do m3=1,3
                        fn(0,konl(m1))=fn(0,konl(m1))-
     &                       c1*qflux(m3)*shp(m3,m1)
                     enddo
                  enddo
               endif
            endif
         enddo
!
!        q contains the contributions to the nodal force in the nodes
!        belonging to the element at stake from other elements (elements
!        already treated). These contributions have to be
!        subtracted to get the contributions attributable to the element
!        at stake only
!
         if(calcul_qa.eq.1) then
            do m1=1,nope
               qa(2)=qa(2)+dabs(fn(0,konl(m1))-q(0,m1))
            enddo
            nal=nal+nope
         endif
      enddo
!
      return
      end
