!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine e_c3d_v1rhs(co,nk,konl,lakonl,p1,p2,omx,bodyfx,
     &     nbody,ff,nelem,nmethod,rhcon,nrhcon,ielmat,ntmat_,vold,vcon,
     &     dtimef,mi,ttime,time,istep,shcon,nshcon,
     &     iturbulent,nelemface,sideface,nface,compressible,
     &     ipvar,var,ipvarf,varf,ithermal,ipface,nelemload,
     &     sideload,xload,nload,ifreesurface,depth,dgravity,cocon,
     &     ncocon,ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iinc,
     &     theta1,bb,physcon,reltimef,xloadold)
!     
!     computation of the velocity element matrix and rhs for the element with
!     element with the topology in konl: step 1 (correction *)
!     
!     ff: rhs 
!     
      implicit none
!     
      character*1 sideface(*)
      character*8 lakonl
      character*20 sideload(*)
!     
      integer konl(8),ifaceq(8,6),nk,nbody,nelem,ithermal(*),mi(*),
     &     i,j,k,i1,i2,j1,nmethod,ii,jj,id,ipointer,idf,
     &     ig,kk,nrhcon(*),ielmat(mi(3),*),nshcon(*),ntmat_,nope,
     &     nopes,imat,iturbulent,compressible,ipface(*),nelemload(2,*),
     &     mint2d,mint3d,ifacet(6,4),ifacew(8,5),istep,nload,
     &     k1,nelemface(*),nface,ipvar(*),index,ipvarf(*),igl,
     &     ifreesurface,ncocon(2,*),ipompc(*),nodempc(3,*),nmpc,
     &     ikmpc(*),ilmpc(*),iinc,iscale,jltyp,nfield,iemchange,
     &     iflux,node
!     
      real*8 co(3,*),shp(4,8),dvi,p1(3),p2(3),x2d3,dep,dvel,
     &     bodyfx(3),ff(0:mi(2),8),bf(3),q(3),c1,c2,xsjmod,dtimef2,fric,
     &     rhcon(0:1,ntmat_,*),vel(3),div,shcon(0:3,ntmat_,*),cp,
     &     voldl(0:mi(2),8),xsj2(3),shp2(7,8),omcor,shps(8),ctuf,
     &     vold(0:mi(2),*),om,omx,const,xsj,temp,tt(3),areaj,enthalpy,
     &     vcon(nk,0:mi(2)),vconl(0:mi(2),8),rho,shpv(8),t(3,3),
     &     cvel(3),vkl(3,3),corio(3),xkin,umttot,xload(2,*),f1,f1m,
     &     xtuf,vort,y,f2,unt,a1,arg2,arg1,gamm,
     &     var(*),varf(*),tu,depth(*),dgravity,gamm1,
     &     beta,beta1,beta2,betas,bfv,c3,c4,cdktuf,ckin,cond,gamm2,
     &     cocon(0:6,ntmat_,*),coefmpc(*),press,skin,skin1,skin2,stuf,
     &     stuf1,stuf2,theta1,tuk,turbprandl,tut,umsk,umst,umt,xkappa,
     &     xtu,tvar(2),rhovel(3),dtem(3),dpress(3),aux(3),coords(3),
     &     bb(3,8),tv(3),pgauss(3),physcon(*),dxkin(3),dxtuf(3),
     &     dxsj2,field,reltimef,sinktemp,tvn,xloadold(2,*),tvnk,tvnt,
     &     xsjmodk,xsjmodt,tvk(3),tvt(3)
!     
      real*8 dtimef,ttime,time
!     
      ifaceq=reshape((/4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/),(/8,6/))
      ifacet=reshape((/1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/),(/6,4/))
      ifacew=reshape((/1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/),(/8,5/))
!     
      tvar(1)=time
      tvar(2)=ttime+dtimef
!     
!     parameter to express 2.d0/3.d0
!     
      x2d3=2.d0/3.d0
      dtimef2=dtimef/2.d0
!     
!     turbulence constants (SST: iturbulent=4)
!
      if(iturbulent.gt.0) then
        a1=0.31d0
        if(iturbulent.eq.4) then
          skin1=0.85d0
        else
          skin1=0.5d0
        endif
        skin2=1.d0
        stuf1=0.5d0
        stuf2=0.856d0
        beta1=0.075d0
        beta2=0.0828d0
        betas=0.09d0
        xkappa= 0.41d0
!     
        gamm1=beta1/betas-stuf1*xkappa*xkappa/dsqrt(betas)
        gamm2=beta2/betas-stuf2*xkappa*xkappa/dsqrt(betas)
!     
        xtu=10.d0*physcon(5)/physcon(8)
        xtu=xtu*xtu
        c3=betas*xtu*10.d0**(-3.5d0)
      endif
!     
      imat=ielmat(1,nelem)
!     
      if(lakonl(4:4).eq.'4') then
        nope=4
        mint3d=1
      elseif(lakonl(4:4).eq.'6') then
        nope=6
        mint3d=2
      elseif(lakonl(4:5).eq.'8R') then
        nope=8
        mint3d=1
      elseif(lakonl(4:4).eq.'8') then
        nope=8
        mint3d=8
      endif
!     
!     initialisation for distributed forces
!     
      do i1=1,nope
        do i2=0,mi(2)
          ff(i2,i1)=0.d0
        enddo
      enddo
      do i1=1,nope
        do i2=1,3
          bb(i2,i1)=0.d0
        enddo
      enddo
!     
!     temperature, velocity, conservative variables
!     (rho*velocity and rho) and if turbulent 
!     rho*turbulence variables
!     
      do i1=1,nope
        do i2=0,mi(2)
          voldl(i2,i1)=vold(i2,konl(i1))
        enddo
        do i2=0,mi(2)
          vconl(i2,i1)=vcon(konl(i1),i2)
        enddo
c        if(nelem.eq.1) then
c          write(*,*) 'voldl(5..6) ',voldl(5,i1),voldl(6,i1)
c          write(*,*) 'vconl(5..6) ',vconl(5,i1),vconl(6,i1)
c        endif
      enddo
!     
!     computation of the matrix: loop over the Gauss points
!     
      index=ipvar(nelem)
      do kk=1,mint3d
!     
!     copying the shape functions, their derivatives and the
!     Jacobian determinant from field var
!     
        do jj=1,nope
          do ii=1,4
            index=index+1
            shp(ii,jj)=var(index)
          enddo
        enddo
        index=index+1
        xsj=var(index)
        index=index+1
        y=var(index)
!     
        xsjmod=dtimef*xsj
!     
!     the temperature temp
!     the velocity vel(*)
!     rho times the velocity rhovel(*)
!     the temperature gradient dtem(*)        
!     the velocity gradient vkl(*,*)
!     the pressure gradient dpress(*)
!     
        temp=0.d0
        do j1=1,3
          vel(j1)=0.d0
          rhovel(j1)=0.d0
          dtem(j1)=0.d0
          do k1=1,3
            vkl(j1,k1)=0.d0
          enddo
          dpress(j1)=0.d0
        enddo
!     
        do i1=1,nope
          temp=temp+shp(4,i1)*voldl(0,i1)
          do j1=1,3
            vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
            rhovel(j1)=rhovel(j1)+shp(4,i1)*vconl(j1,i1)
            dtem(j1)=dtem(j1)+shp(j1,i1)*voldl(0,i1)
            do k1=1,3
              vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
            enddo
            dpress(j1)=dpress(j1)+shp(j1,i1)*voldl(4,i1)
          enddo
        enddo
!     
!     the divergence of the velocity div 
!
        if(compressible.eq.1) then
          div=vkl(1,1)+vkl(2,2)+vkl(3,3)
        else
          div=0.d0
        endif
!     
!     the divergence of the shape function times the velocity shpv(*)
!     the convective enthalpy        
!     the convective velocity cvel(*)
!     
        shpv(1)=0.d0
        enthalpy=0.d0
        do j1=1,3
          cvel(j1)=0.d0
        enddo
!     
        do i1=1,nope
          shpv(i1)=shp(1,i1)*vel(1)+shp(2,i1)*vel(2)+
     &         shp(3,i1)*vel(3)+shp(4,i1)*div
          enthalpy=enthalpy+shpv(i1)*(vconl(0,i1)+voldl(4,i1))
          do j1=1,3
            cvel(j1)=cvel(j1)+shpv(i1)*vconl(j1,i1)
          enddo
        enddo
!     
!     creating auxiliary variable shps
!     
        do i1=1,nope
          shps(i1)=xsjmod*(shp(4,i1)+dtimef2*shpv(i1))
        enddo
!     
!     material data (dynamic viscosity)
!     
        call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,dvi)
!     
!     determining the dissipative stress 
!     
        do i1=1,3
          do j1=i1,3
            t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
          enddo
          if(compressible.eq.1) t(i1,i1)=t(i1,i1)-x2d3*div
        enddo
!     
!     calculation of the density in case of body forces and/or
!     turbulence
!     
        if((nbody.ne.0).or.(iturbulent.ne.0)) then
          if(compressible.eq.1) then
!     
!     gas
!     
            rho=0.d0
            do i1=1,nope
              rho=rho+shp(4,i1)*vconl(4,i1)
            enddo
          else
!     
!     liquid
!     
            call materialdata_rho(rhcon,nrhcon,imat,rho,
     &           temp,ntmat_,ithermal)
          endif
        endif
!
!     calculation of the turbulent kinetic energy, turbulence
!     frequency and their spatial derivatives for gases and liquids
!
        if(iturbulent.gt.0) then
          xkin=0.d0
          xtuf=0.d0
          do i1=1,nope
            xkin=xkin+shp(4,i1)*voldl(5,i1)
            xtuf=xtuf+shp(4,i1)*voldl(6,i1)
          enddo
!     
!     adding the turbulent stress
!     
!     factor F2
!     
          c1=dsqrt(xkin)/(0.09d0*xtuf*y)
          c2=500.d0*dvi/(y*y*xtuf*rho)
!     
!     kinematic turbulent viscosity
!     
          if(iturbulent.eq.4) then
!     
!     vorticity
!     
            vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &           (vkl(1,3)-vkl(3,1))**2+
     &           (vkl(2,1)-vkl(1,2))**2)
            arg2=max(2.d0*c1,c2)
            f2=dtanh(arg2*arg2)
            unt=a1*xkin/max(a1*xtuf,vort*f2)
          else
            unt=xkin/xtuf
          endif
!     
!     calculating the production (anisotropic part of
!     the turbulent stress is, apart from the dynamic
!     viscosity, identical to the viscous stress)
!     
          tu=t(1,1)*vkl(1,1)+t(1,2)*t(1,2)+t(1,3)*t(1,3)+
     &         t(2,2)*vkl(2,2)+t(2,3)*t(2,3)+t(3,3)*vkl(3,3)
!     
!     correction for compressible fluids
!     
          if(compressible.eq.1) then
            tu=tu-x2d3*xkin*div/unt
          endif
!     
!     calculating the turbulent stress
!     
c          umttot=dvi
          umttot=dvi+unt*rho
          do i1=1,3
            do j1=i1,3
              t(i1,j1)=umttot*t(i1,j1)
            enddo
            t(i1,i1)=t(i1,i1)-x2d3*rho*xkin
          enddo
        else
!     
          do i1=1,3
            do j1=i1,3
              t(i1,j1)=dvi*t(i1,j1)
            enddo
          enddo
        endif
!     
!     1a. rhs of the first part of the momentum equation
!     
        do jj=1,nope
!     
!     convective + diffusive
!     
          ff(1,jj)=ff(1,jj)-cvel(1)*shps(jj)-xsjmod*
     &         (shp(1,jj)*t(1,1)+shp(2,jj)*t(1,2)+shp(3,jj)*t(1,3))
          ff(2,jj)=ff(2,jj)-cvel(2)*shps(jj)-xsjmod*
     &         (shp(1,jj)*t(1,2)+shp(2,jj)*t(2,2)+shp(3,jj)*t(2,3))
          ff(3,jj)=ff(3,jj)-cvel(3)*shps(jj)-xsjmod*
     &         (shp(1,jj)*t(1,3)+shp(2,jj)*t(2,3)+shp(3,jj)*t(3,3))
        enddo
!     
!     computation of contribution due to body forces
!     
        if(nbody.ne.0) then
!     
!     initialisation for the body forces
!     
          om=omx*rho
          omcor=2.d0*rho*dsqrt(omx)
!     
          if(om.gt.0.d0) then
            do i1=1,3
!     
!     computation of the global coordinates of the gauss
!     point
!     
              q(i1)=0.d0
              do j1=1,nope
                q(i1)=q(i1)+shp(4,j1)*co(i1,konl(j1))
              enddo
!     
              q(i1)=q(i1)-p1(i1)
            enddo
            const=q(1)*p2(1)+q(2)*p2(2)+q(3)*p2(3)
!     
!     Coriolis forces
!     
            omcor=2.d0*rho*dsqrt(omx)
            corio(1)=vel(2)*p2(3)-vel(3)*p2(2)
            corio(2)=vel(3)*p2(1)-vel(1)*p2(3)
            corio(3)=vel(1)*p2(2)-vel(2)*p2(1)
          endif
!     
          if(ifreesurface.eq.0) then
            do ii=1,3
              bf(ii)=bodyfx(ii)*rho
            enddo
!     
!     inclusion of the centrifugal force into the body force
!     
            if(om.gt.0.d0) then
              do i1=1,3
                bf(i1)=bf(i1)+(q(i1)-const*p2(i1))*om+
     &               corio(i1)*omcor
              enddo
            endif
          else
!     
!     shallow water calculation
!     effect of varying depth;
!     the effect of the centrifugal force on dgravity is neglected            
!     
            dep=0.d0
            do j1=1,3
              bf(j1)=0.d0
            enddo
!     
            do i1=1,nope
              dep=dep+shp(4,i1)*depth(konl(i1))
              do j1=1,3
                bf(j1)=bf(j1)+shp(j1,i1)*depth(konl(i1))
              enddo
            enddo
            do j1=1,3
              bf(j1)=bf(j1)*(rho-dep)*dgravity
            enddo
!     
            if(om.gt.0.d0) then
              do i1=1,2
                bf(i1)=bf(i1)+(q(i1)-const*p2(i1))*om+
     &               corio(i1)*omcor
              enddo
            endif
!     
!     bottom friction
!     
c     fric=0.02d0
            fric=0.01d0
            dvel=dsqrt(vel(1)*vel(1)+vel(2)*vel(2)+vel(3)*vel(3))
            do j1=1,3
              bf(j1)=bf(j1)-fric*dvel*vel(j1)/8.d0
            enddo
          endif
!     
!     storing the body force
!     
          bfv=bf(1)*vel(1)+bf(2)*vel(2)+bf(3)*vel(3)
!     
!     1b. rhs of the first part of the momentum equation:
!         body force contribution
!     
          do jj=1,nope
            ff(1,jj)=ff(1,jj)+bf(1)*shps(jj)
            ff(2,jj)=ff(2,jj)+bf(2)*shps(jj)
            ff(3,jj)=ff(3,jj)+bf(3)*shps(jj)
          enddo
        endif
!
!       2. rhs of the mass equation
!
        do j1=1,3
          aux(j1)=xsjmod*(rhovel(j1)-dtimef*theta1*dpress(j1))
        enddo
!        
        do jj=1,nope
          ff(4,jj)=ff(4,jj)+
     &         shp(1,jj)*aux(1)+shp(2,jj)*aux(2)+shp(3,jj)*aux(3)
        enddo
!     
!     3. rhs of the second part of the momentum equation:
!     
        if(compressible.eq.1) then
!
!         explicit compressible
!
          do jj=1,nope
            bb(1,jj)=bb(1,jj)-dpress(1)*shps(jj)
            bb(2,jj)=bb(2,jj)-dpress(2)*shps(jj)
            bb(3,jj)=bb(3,jj)-dpress(3)*shps(jj)
          enddo
        else
!
!         implicit incompressible
!
          do jj=1,nope
            bb(1,jj)=bb(1,jj)-xsjmod*shp(4,jj)*dpress(1)
            bb(2,jj)=bb(2,jj)-xsjmod*shp(4,jj)*dpress(2)
            bb(3,jj)=bb(3,jj)-xsjmod*shp(4,jj)*dpress(3)
          enddo
        endif
!     
!     4. rhs of the energy equation:
!     
        if(ithermal(1).gt.0) then
!     
!     viscous conductivity
!     
          call materialdata_cond(imat,ntmat_,temp,cocon,ncocon,cond)
!     
!     adding the turbulent conductivity
!     
          if(iturbulent.gt.0) then
            call materialdata_cp(imat,ntmat_,temp,shcon,nshcon,cp)
            turbprandl=0.9d0
            cond=cond+cp*rho*unt/turbprandl
          endif
!     
!     calculating the total dissipative stress x velocity
!     (viscous + turbulent)
!     
          tv(1)=t(1,1)*vel(1)+t(1,2)*vel(2)+t(1,3)*vel(3)
          tv(2)=t(1,2)*vel(1)+t(2,2)*vel(2)+t(2,3)*vel(3)
          tv(3)=t(1,3)*vel(1)+t(2,3)*vel(2)+t(3,3)*vel(3)
!     
!     determining stress x velocity + conductivity x
!     temperature gradient
!     
          do i1=1,3
            tv(i1)=tv(i1)+cond*dtem(i1)
          enddo
!     
!     determination of the rhs of the energy equations
!     
          do jj=1,nope
            ff(0,jj)=ff(0,jj)-shps(jj)*enthalpy-xsjmod*
     &           (shp(1,jj)*tv(1)+shp(2,jj)*tv(2)+shp(3,jj)*tv(3))
          enddo
!     
!     computation of contribution due to body forces
!     
          if(nbody.ne.0) then
            do jj=1,nope
              ff(0,jj)=ff(0,jj)+shps(jj)*bfv
            enddo
          endif
!     
!     distributed heat flux
!     
          if(nload.gt.0) then
            call nident2(nelemload,nelem,nload,id)
            areaj=xsj
            do
              if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
              if(sideload(id)(1:2).ne.'BF') then
                id=id-1
                cycle
              endif
              if(sideload(id)(3:4).eq.'NU') then
                do j=1,3
                  pgauss(j)=0.d0
                  do i1=1,nope
                    pgauss(j)=pgauss(j)+
     &                   shp(4,i1)*co(j,konl(i1))
                  enddo
                enddo
                jltyp=1
                iscale=1
                call dflux(xload(1,id),temp,istep,iinc,tvar,
     &               nelem,kk,pgauss,jltyp,temp,press,sideload(id),
     &               areaj,vold,co,lakonl,konl,ipompc,nodempc,coefmpc,
     &               nmpc,ikmpc,ilmpc,iscale,mi)
              endif
              do jj=1,nope
                ff(0,jj)=ff(0,jj)+shps(jj)*xload(1,id)
              enddo
              exit
            enddo
          endif
        endif
!     
!     5. rhs of the turbulence equations:
!     
        if(iturbulent.gt.0) then
!     
!     convective turbulent kinetic energy: ckin
!     convective turbulence frequency: ctuf
!     
          ckin=0.d0
          ctuf=0.d0
          do i1=1,nope
            ckin=ckin+shpv(i1)*vconl(5,i1)
            ctuf=ctuf+shpv(i1)*vconl(6,i1)
          enddo
c        if(nelem.eq.1) then
c          write(*,*) 'ckin ',ckin
c          write(*,*) 'ctuf ',ctuf
c        endif
!     
!     gradient of k and omega
!     
          do j1=1,3
            dxkin(j1)=0.d0
            dxtuf(j1)=0.d0
          enddo
          do i1=1,nope
            do j1=1,3
              dxkin(j1)=dxkin(j1)+shp(j1,i1)*voldl(5,i1)
              dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*voldl(6,i1)
            enddo
          enddo
c        if(nelem.eq.1) then
c          write(*,*) 'dxkin ',(dxkin(j1),j1=1,3)
c          write(*,*) 'dxtuf ',(dxtuf(j1),j1=1,3)
c        endif
!     
!     auxiliary variable
!     
          c4=2.d0*rho*stuf2*
     &         (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+dxkin(3)*dxtuf(3))/
     &         xtuf
!     
!     dynamic turbulent viscosity
!     
          umt=unt*rho
!     
!     factor F1
!     
          if(iturbulent.eq.1) then
!     
!     k-epsilon model
!     
            f1=0.d0
          elseif(iturbulent.eq.2) then
!     
!     k-omega model
!     
            f1=1.d0
          else
!     
!     BSL/SST model
!     
            cdktuf=max(c4,1.d-20)
            arg1=min(max(c1,c2),4.d0*rho*stuf2*xkin/(cdktuf*y*y))
            f1=dtanh(arg1**4.d0)
          endif
          f1m=1.d0-f1
!     
!     interpolation of the constants
!     
          skin=f1*skin1+f1m*skin2
          stuf=f1*stuf1+f1m*stuf2
          beta=f1*beta1+f1m*beta2
          gamm=f1*gamm1+f1m*gamm2
!     
!     source terms: productivity - dissipation
!     
          umsk=dvi+skin*umt
          umst=dvi+stuf*umt
!     
!     production limiter active: P=unt*tu<=20*betas*k*omega
!     Menter, F.R., "Zonal Two Equation k-omega Turbulence Models for
!     Aerodynamic Flows," AIAA Paper 93-2906, July 1993.
!     
          tuk=rho*(unt*tu-betas*xtuf*xkin)
          tut=rho*(gamm*tu-beta*xtuf*xtuf)+f1m*c4
!     
!     add controlled decay
!     Spalart, P.R. and Rumsey, C.L., "Effective Inflow Conditions for
!     Turbulence Models in Aerodynamic Calculations," AIAA Journal,
!     Vol. 45, No. 10, 2007,pp.2544-2553.
!     
          tuk=tuk+c3*dvi
          tut=tut+beta*xtu*rho
!     
          do i1=1,3
            dxkin(i1)=dxkin(i1)*umsk
            dxtuf(i1)=dxtuf(i1)*umst
          enddo
!     
!     determination of rhs
!     
          do jj=1,nope
!     
            ff(5,jj)=ff(5,jj)-shps(jj)*(ckin-tuk)-xsjmod*
     &           (shp(1,jj)*dxkin(1)+shp(2,jj)*dxkin(2)
     &           +shp(3,jj)*dxkin(3))
            ff(6,jj)=ff(6,jj)-shps(jj)*(ctuf-tut)-xsjmod*
     &           (shp(1,jj)*dxtuf(1)+shp(2,jj)*dxtuf(2)
     &           +shp(3,jj)*dxtuf(3))
          enddo
        endif
!        
      enddo
!
!     area integrals
!     
      if(nface.ne.0) then
        index=ipvarf(nelem)
c        write(*,*) 'e_c3d_v1rhs nelem ipvarf(nelem)',nelem,index
!     
!     external boundaries
!     
        nopes=0
        idf=ipface(nelem)
        do
          if((idf.eq.0).or.(nelemface(idf).ne.nelem)) exit
          ig=ichar(sideface(idf)(1:1))-48
!     
!     check for distributed flux
!     an adiabatic face must be declared as a face with
!     distributed flux zero!
!     
          iflux=0
          call nident2(nelemload,nelem,nload,id)
          do
            if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
            if((sideload(id)(1:1).ne.'F').and.
     &           (sideload(id)(1:1).ne.'R').and.
     &           (sideload(id)(1:1).ne.'S')) then
              id=id-1
              cycle
            endif
            igl=ichar(sideload(id)(2:2))-48
            if(igl.ne.ig) then
              id=id-1
              cycle
            endif
            iflux=1
            exit
          enddo
!     
          if(nopes.eq.0) then
            if(lakonl(4:4).eq.'4') then
              nopes=3
              mint2d=1
            elseif(lakonl(4:4).eq.'6') then
              mint2d=1
            elseif(lakonl(4:5).eq.'8R') then
              nopes=4
              mint2d=1
            elseif(lakonl(4:4).eq.'8') then
              nopes=4
              mint2d=4
            endif
          endif
!     
          if(lakonl(4:4).eq.'6') then
            if(ig.le.2) then
              nopes=3
            else
              nopes=4
            endif
          endif
!
c          write(*,*) 'e_c3d_v1rhs ',index,4*nope+nopes+4
          do i=1,mint2d
!     
!     facial shape functions
!     local surface normal
!     
            do i1=1,nopes
              index=index+1
              shp2(4,i1)=varf(index)
            enddo
            do i1=1,3
              index=index+1
              xsj2(i1)=varf(index)
            enddo
!     
!     derivative of the volumetric shape functions
!     needed for the temperature, velocity gradients and
!     gradients of k and omega (turbulence)
!     
            do i1=1,nope
              do j1=1,4
                index=index+1
                shp(j1,i1)=varf(index)
              enddo
            enddo
            index=index+1
            y=varf(index)
!     
!     calculating of
!     the temperature temp
!     the velocity vel(*)
!     rho times the velocity rhovel(*)
!     the velocity gradient vkl
!     
            temp=0.d0
            do j1=1,3
              vel(j1)=0.d0
              rhovel(j1)=0.d0
              do k1=1,3
                vkl(j1,k1)=0.d0
              enddo
            enddo
!
            do i1=1,nope
              temp=temp+shp(4,i1)*voldl(0,i1)
              do j1=1,3
                vel(j1)=vel(j1)+shp(4,i1)*voldl(j1,i1)
c                write(*,*) 'e_c3d_v1rhs ',shp(4,i1)
                rhovel(j1)=rhovel(j1)+shp(4,i1)*vconl(j1,i1)
                do k1=1,3
                  vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
                enddo
              enddo
            enddo
!     
            if(iflux.eq.0) then
!     
!     calculating of the temperature gradient dtem
!     in the integration point
!
              do j1=1,3
                dtem(j1)=0.d0
              enddo
              do i1=1,nope
                do j1=1,3
                  dtem(j1)=dtem(j1)+shp(j1,i1)*voldl(0,i1)
                enddo
              enddo
            endif
!     
            if(compressible.eq.1) then
              div=vkl(1,1)+vkl(2,2)+vkl(3,3)
            else
              div=0.d0
            endif
!     
!     material data (dynamic viscosity)
!
            call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,
     &           dvi)
!     
!     determining the dissipative stress 
!     
            do i1=1,3
              do j1=i1,3
                t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
              enddo
              if(compressible.eq.1) t(i1,i1)=t(i1,i1)-x2d3*div
            enddo
!     
!     calculation of the density for gases
!     
!     calculation of the turbulent kinetic energy, turbulence
!     frequency and their spatial derivatives for gases and liquids
!     
            if(iturbulent.gt.0) then
              if(compressible.eq.1) then
!     
!     gas
!     
                rho=0.d0
                do i1=1,nope
                  rho=rho+shp(4,i1)*vconl(4,i1)
                enddo
              else
!     
!     liquid
!     
                call materialdata_rho(rhcon,nrhcon,imat,rho,
     &               temp,ntmat_,ithermal)
!
!     calculation of k, omega an y
!
              endif
              xkin=0.d0
              xtuf=0.d0
              do i1=1,nope
                xkin=xkin+shp(4,i1)*voldl(5,i1)
                xtuf=xtuf+shp(4,i1)*voldl(6,i1)
              enddo
!     
!     calculation of turbulent auxiliary variables
!     
!     factor F2
!
              if(y.gt.0.d0) then
                c1=dsqrt(xkin)/(0.09d0*xtuf*y)
                c2=500.d0*dvi/(y*y*xtuf*rho)
              endif
!     
!     kinematic and dynamic turbulent viscosity
!     
              if(iturbulent.eq.4) then
!     
!     vorticity
!     
                vort=dsqrt((vkl(3,2)-vkl(2,3))**2+
     &               (vkl(1,3)-vkl(3,1))**2+
     &               (vkl(2,1)-vkl(1,2))**2)
                if(y.gt.0.d0) then
                  arg2=max(2.d0*c1,c2)
                  f2=dtanh(arg2*arg2)
                else
                  f2=1.d0
                endif
                unt=a1*xkin/max(a1*xtuf,vort*f2)
              else
                unt=xkin/xtuf
              endif
!     
              umttot=dvi+unt*rho
              do i1=1,3
                do j1=i1,3
                  t(i1,j1)=umttot*t(i1,j1)
                enddo
                t(i1,i1)=t(i1,i1)-x2d3*rho*xkin
              enddo
            else
!     
              do i1=1,3
                do j1=i1,3
                  t(i1,j1)=dvi*t(i1,j1)
                enddo
              enddo
            endif
!
!     stress vector
!
            tt(1)=(t(1,1)*xsj2(1)+t(1,2)*xsj2(2)+t(1,3)*xsj2(3))*
     &           dtimef
            tt(2)=(t(1,2)*xsj2(1)+t(2,2)*xsj2(2)+t(2,3)*xsj2(3))*
     &           dtimef
            tt(3)=(t(1,3)*xsj2(1)+t(2,3)*xsj2(2)+t(3,3)*xsj2(3))*
     &           dtimef
!     
!      stress x velocity
!     
            tv(1)=t(1,1)*vel(1)+t(1,2)*vel(2)+t(1,3)*vel(3)
            tv(2)=t(1,2)*vel(1)+t(2,2)*vel(2)+t(2,3)*vel(3)
            tv(3)=t(1,3)*vel(1)+t(2,3)*vel(2)+t(3,3)*vel(3)
!
!     adding conductivity in case the flux is not given by a
!     *DFLUX, *FILM or *RADIATE card
!     
            if(iflux.eq.0) then
              do i1=1,3
                tv(i1)=tv(i1)+cond*dtem(i1)
              enddo
            endif
!     
            tvn=tv(1)*xsj2(1)+tv(2)*xsj2(2)+tv(3)*xsj2(3)
!
!     modifying tvn in case of a *DFLUX, *FILM or *RADIATE
!     card
!
            if(iflux.eq.1) then
              dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &             xsj2(3)*xsj2(3))
              areaj=dxsj2
              sinktemp=xload(2,id)
!     
!     for nonuniform load: determine the coordinates of the
!     point (transferred into the user subroutine)
!     
              if((sideload(id)(3:4).eq.'NU').or.
     &             (sideload(id)(5:6).eq.'NU')) then
                if(nope.eq.8) then
                  do k=1,3
                    coords(k)=0.d0
                    do j=1,nopes
                      coords(k)=coords(k)+
     &                     co(k,konl(ifaceq(j,ig)))*shp2(4,j)
                    enddo
                  enddo
                elseif(nope.eq.4) then
                  do k=1,3
                    coords(k)=0.d0
                    do j=1,nopes
                      coords(k)=coords(k)+
     &                     co(k,konl(ifacet(j,ig)))*shp2(4,j)
                    enddo
                  enddo
                else
                  do k=1,3
                    coords(k)=0.d0
                    do j=1,nopes
                      coords(k)=coords(k)+
     &                     co(k,konl(ifacew(j,ig)))*shp2(4,j)
                    enddo
                  enddo
                endif
                jltyp=ichar(sideload(id)(2:2))-48
                jltyp=jltyp+10
                if(sideload(id)(1:1).eq.'S') then
                  iscale=1
                  call dflux(xload(1,id),temp,istep,iinc,tvar,
     &                 nelem,i,coords,jltyp,temp,press,
     &                 sideload(id),areaj,vold,co,lakonl,konl,
     &                 ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,
     &                 iscale,mi)
                  if((nmethod.eq.1).and.(iscale.ne.0))
     &                 xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltimef
                elseif(sideload(id)(1:1).eq.'F') then
                  call film(xload(1,id),sinktemp,temp,istep,
     &                 iinc,tvar,nelem,i,coords,jltyp,field,
     &                 nfield,sideload(id),node,areaj,vold,mi)
                  if(nmethod.eq.1) xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltimef
                elseif(sideload(id)(1:1).eq.'R') then
                  call radiate(xload(1,id),xload(2,id),temp,istep,
     &                 iinc,tvar,nelem,i,coords,jltyp,field,
     &                 nfield,sideload(id),node,areaj,vold,mi,
     &                 iemchange)
                  if(nmethod.eq.1) xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltimef
                endif
              endif
!     
              if(sideload(id)(1:1).eq.'S') then
!     
!     flux INTO the face is positive (input deck convention)
!     this is different from the convention in the theory
!     
                tvn=tvn+xload(1,id)*dxsj2
              elseif(sideload(id)(1:1).eq.'F') then
                tvn=tvn-xload(1,id)*(temp-sinktemp)*dxsj2
              elseif(sideload(id)(1:1).eq.'R') then
                tvn=tvn-physcon(2)*
     &               xload(1,id)*((temp-physcon(1))**4-
     &               (xload(2,id)-physcon(1))**4)*dxsj2
              endif
            endif
!     
            xsjmod=tvn*dtimef
!
            if(iturbulent.gt.0) then
!     
!     calculation of the spatial derivatives of the turbulent kinetic energy
!     and the turbulence frequency for gases and liquids
!     
              do j1=1,3
                dxkin(j1)=0.d0
                dxtuf(j1)=0.d0
              enddo
              do i1=1,nope
                do j1=1,3
                  dxkin(j1)=dxkin(j1)+shp(j1,i1)*voldl(5,i1)
                  dxtuf(j1)=dxtuf(j1)+shp(j1,i1)*voldl(6,i1)
                enddo
              enddo
!     
!     auxiliary variable
!     
              c4=2.d0*rho*stuf2*
     &             (dxkin(1)*dxtuf(1)+dxkin(2)*dxtuf(2)+
     &             dxkin(3)*dxtuf(3))/xtuf
!     
!     dynamic turbulent viscosity
!     
              umt=unt*rho
!     
!     factor F1
!                
              if(iturbulent.eq.1) then
!     
!     k-epsilon model
!     
                f1=0.d0
              elseif(iturbulent.eq.2) then
!     
!     k-omega model
!     
                f1=1.d0
              else
!     
!     BSL/SST model
!     
                if(y.gt.0.d0) then
!     
!     finite distance from wall
!     
                  cdktuf=max(c4,1.d-20)
                  arg1=
     &                 min(max(c1,c2),4.d0*rho*stuf2*xkin/(cdktuf*y*y))
                  f1=dtanh(arg1**4.d0)
                else
!     
!     wall
!     
                  f1=1.d0
                endif
              endif
              f1m=1.d0-f1
!     
!     interpolation of the constants
!     
              skin=f1*skin1+f1m*skin2
              stuf=f1*stuf1+f1m*stuf2
!     
!     auxiliary quantities
!     
              umsk=dvi+skin*umt
              umst=dvi+stuf*umt
!     
!     determining the stress and and stress x velocity + conductivity x
!     temperature gradient
!     
              do i1=1,3
                tvk(i1)=umsk*dxkin(i1)
                tvt(i1)=umsk*dxtuf(i1)
              enddo
!     
              tvnk=tvk(1)*xsj2(1)+tvk(2)*xsj2(2)+tvk(3)*xsj2(3)
              tvnt=tvt(1)*xsj2(1)+tvt(2)*xsj2(2)+tvt(3)*xsj2(3)
!     
              xsjmodk=tvnk*dtimef
              xsjmodt=tvnt*dtimef
            endif
!
            do k=1,nopes
              if(nope.eq.8) then
                ipointer=ifaceq(k,ig)
              elseif(nope.eq.4) then
                ipointer=ifacet(k,ig)
              else
                ipointer=ifacew(k,ig)
              endif
!     
!     1a. rhs of the first part of the momentum equation
!     
              ff(1,ipointer)=ff(1,ipointer)+shp2(4,k)*tt(1)
              ff(2,ipointer)=ff(2,ipointer)+shp2(4,k)*tt(2)
              ff(3,ipointer)=ff(3,ipointer)+shp2(4,k)*tt(3)
!
!     2. rhs of the mass equation
!
              ff(4,ipointer)=ff(4,ipointer)-shp2(4,k)*
     &             (rhovel(1)*xsj2(1)+rhovel(2)*xsj2(2)+
     &             rhovel(3)*xsj2(3))*dtimef
!     
!     4. rhs of the energy equation:
!     
              if(ithermal(1).gt.0) then
                ff(0,ipointer)=ff(0,ipointer)+shp2(4,k)*xsjmod
              endif
!     
!     5. rhs of the turbulence equations:
!
              if(iturbulent.gt.0) then
                ff(5,ipointer)=ff(5,ipointer)+shp2(4,k)*xsjmodk
                ff(6,ipointer)=ff(6,ipointer)+shp2(4,k)*xsjmodt
              endif
            enddo
          enddo
!
          idf=idf-1
        enddo
      endif
!
c      if(nelem.eq.1) then
c      write(*,*) 'e_c3d_v1rhs'
c      do k=1,8
c        write(*,*) nelem,k,(ff(j,k),j=5,6)
c      enddo
c        endif
      return
      end
