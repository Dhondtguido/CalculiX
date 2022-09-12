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
      subroutine printoutface(co,ntmat_,v,
     &  cocon,ncocon,istartset,iendset,ipkon,lakon,kon,
     &  ialset,prset,ttime,nset,set,nprint,prlab,ielmat,mi,
     &  time,stn,iperturb)
!
!     calculation and printout of the heat flux, forces and/or
!     moments on a surface
!
      implicit none
!
      character*8 lakonl,lakon(*)
      character*6 prlab(*)
      character*80 faset
      character*81 set(*),prset(*)
!
      integer konl(20),ifaceq(8,6),nelem,ii,nprint,i,j,i1,i2,j1,
     &  ncocon(2,*),k1,jj,ig,ntmat_,nope,nopes,imat,iperturb(*),
     &  mint2d,ifacet(6,4),ifacew(8,5),iflag,indexe,jface,istartset(*),
     &  iendset(*),ipkon(*),kon(*),iset,ialset(*),nset,ipos,id,
     &  mi(*),ielmat(mi(3),*)
!
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),f(0:3),time,
     &  vkl(0:3,3),t(3,3),xtorque,bendmom,xnormforc,shearforc,
     &  vl(0:mi(2),20),cocon(0:6,ntmat_,*),xl2(3,8),xsj2(3),
     &  shp2(7,8),v(0:mi(2),*),xi,et,xsj,temp,xi3d,et3d,ze3d,weight,
     &  xlocal20(3,9,6),xlocal4(3,1,4),xlocal10(3,3,4),xlocal6(3,1,5),
     &  xlocal15(3,4,5),xlocal8(3,4,6),xlocal8r(3,1,6),ttime,
     &  tf(0:3),dd,coords(3),cond,stn(6,*),xm(3),df(3),cg(3),
     &  area,xn(3),xmcg(3)
!
      include "gauss.f"
      include "xlocal.f"
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data iflag /3/
!
!     flag to check whether sof and/or som output was already
!     printed
!
      do ii=1,nprint
!
!        check whether there are facial print requests
!
!        FLUX: heat flux
!        SOF: forces and moments  on a section 
!              (= an internal or external surface)
!
!     DRAG removed on 13 Dec 2020: DRAG only accessible for FEM-CBS
!     through printoutfacefem.f
!
         if((prlab(ii)(1:4).eq.'FLUX').or.
     &      (prlab(ii)(1:3).eq.'SOF'))      
     &      then
!
            ipos=index(prset(ii),' ')
            do i=1,80
               faset(i:i)=' '
            enddo
            faset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
!     printing the header
!     
            write(5,*)
            if(prlab(ii)(1:4).eq.'FLUX') then
!
!              initialisierung of the flux
!     
               f(0)=0.d0
               write(5,121) faset(1:ipos-2),ttime+time
 121           format(' heat flux at the integration points for set ',
     &                A,' and time ',e14.7)
               write(5,*)
               write(5,125)
 125           format('        el  fa  int heat flux q and             c
     &oordinates')
            else
!
!              initialize total force and total moment
!
               do i=1,3
                  f(i)=0.d0
                  xm(i)=0.d0
                  cg(i)=0.d0
                  xn(i)=0.d0
               enddo
               area=0.d0
            endif
            write(5,*)
!     
!           printing the data
!
            call cident81(set,prset(ii),nset,id)
            iset=nset+1
            if(id.gt.0) then
              if(prset(ii).eq.set(id)) then
                iset=id
              endif
            endif
!
            do jj=istartset(iset),iendset(iset)
!     
               jface=ialset(jj)
!     
               nelem=int(jface/10.d0)
               ig=jface-10*nelem
!     
               indexe=ipkon(nelem)
               lakonl=lakon(nelem)
               imat=ielmat(1,nelem)
!     
               if(lakonl(4:4).eq.'2') then
                  nope=20
                  nopes=8
               elseif(lakonl(4:4).eq.'8') then
                  nope=8
                  nopes=4
               elseif(lakonl(4:5).eq.'10') then
                  nope=10
                  nopes=6
               elseif(lakonl(4:4).eq.'4') then
                  nope=4
                  nopes=3
               elseif(lakonl(4:5).eq.'15') then
                  nope=15
               elseif(lakonl(4:4).eq.'6') then
                  nope=6
               endif
!     
               if(lakonl(4:5).eq.'8R') then
                  mint2d=1
               elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) 
     &             then
                  if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'S').or.
     &                 (lakonl(7:7).eq.'E')) then
                     mint2d=2
                  else
                     mint2d=4
                  endif
               elseif(lakonl(4:4).eq.'2') then
                  mint2d=9
               elseif(lakonl(4:5).eq.'10') then
                  mint2d=3
               elseif(lakonl(4:4).eq.'4') then
                  mint2d=1
               endif
!     
!     local topology
!     
               do i=1,nope
                  konl(i)=kon(indexe+i)
               enddo
!     
!     computation of the coordinates of the local nodes
!     
               do i=1,nope
                  do j=1,3
                     xl(j,i)=co(j,konl(i))
                  enddo
               enddo
!     
!     temperature, displacement (structures) or
!     
               do i1=1,nope
                  do i2=0,mi(2)
                     vl(i2,i1)=v(i2,konl(i1))
                  enddo
               enddo
!     
!     for structural calculations with NLGEOM: adding the displacements
!     to the coordinates
!
               if(iperturb(2).eq.1) then
                 do i=1,nope
                   do j=1,3
                     xl(j,i)=xl(j,i)+vl(j,i)
                   enddo
                 enddo
               endif
!     
!     treatment of wedge faces
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
!     no CFD: deformed structure
!     
               if(iperturb(2).eq.1) then
                 if((nope.eq.20).or.(nope.eq.8)) then
                   do i=1,nopes
                     do j=1,3
                       xl2(j,i)=co(j,konl(ifaceq(i,ig)))
     &                      +v(j,konl(ifaceq(i,ig)))
                     enddo
                   enddo
                 elseif((nope.eq.10).or.(nope.eq.4)) then
                   do i=1,nopes
                     do j=1,3
                       xl2(j,i)=co(j,konl(ifacet(i,ig)))
     &                      +v(j,konl(ifacet(i,ig)))
                     enddo
                   enddo
                 else
                   do i=1,nopes
                     do j=1,3
                       xl2(j,i)=co(j,konl(ifacew(i,ig)))+
     &                      v(j,konl(ifacew(i,ig)))
                     enddo
                   enddo
                 endif
               else
                 if((nope.eq.20).or.(nope.eq.8)) then
                   do i=1,nopes
                     do j=1,3
                       xl2(j,i)=co(j,konl(ifaceq(i,ig)))
                     enddo
                   enddo
                 elseif((nope.eq.10).or.(nope.eq.4)) then
                   do i=1,nopes
                     do j=1,3
                       xl2(j,i)=co(j,konl(ifacet(i,ig)))
                     enddo
                   enddo
                 else
                   do i=1,nopes
                     do j=1,3
                       xl2(j,i)=co(j,konl(ifacew(i,ig)))
                     enddo
                   enddo
                 endif
               endif
!     
               do i=1,mint2d
!     
!     local coordinates of the surface integration
!     point within the surface local coordinate system
!     
                  if((lakonl(4:5).eq.'8R').or.
     &                 ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                     xi=gauss2d1(1,i)
                     et=gauss2d1(2,i)
                     weight=weight2d1(i)
                  elseif((lakonl(4:4).eq.'8').or.
     &                    (lakonl(4:6).eq.'20R').or.
     &                    ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
                     xi=gauss2d2(1,i)
                     et=gauss2d2(2,i)
                     weight=weight2d2(i)
                  elseif(lakonl(4:4).eq.'2') then
                     xi=gauss2d3(1,i)
                     et=gauss2d3(2,i)
                     weight=weight2d3(i)
                  elseif((lakonl(4:5).eq.'10').or.
     &                    ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
                     xi=gauss2d5(1,i)
                     et=gauss2d5(2,i)
                     weight=weight2d5(i)
                  elseif((lakonl(4:4).eq.'4').or.
     &                    ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                     xi=gauss2d4(1,i)
                     et=gauss2d4(2,i)
                     weight=weight2d4(i)
                  endif
!     
!     local surface normal
!     
                  if(nopes.eq.8) then
                     call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                  elseif(nopes.eq.4) then
                     call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
                  elseif(nopes.eq.6) then
                     call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                  else
                     call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
                  endif
!
!                 global coordinates of the integration point
!
                  do j1=1,3
                     coords(j1)=0.d0
                     do i1=1,nopes
                        coords(j1)=coords(j1)+shp2(4,i1)*xl2(j1,i1)
                     enddo
                  enddo
!     
!     local coordinates of the surface integration
!     point within the element local coordinate system
!     
                  if(prlab(ii)(1:4).eq.'FLUX') then
!
!                    deformation gradient is only needed for
!                    FLUX applications
!
                     if(lakonl(4:5).eq.'8R') then
                        xi3d=xlocal8r(1,i,ig)
                        et3d=xlocal8r(2,i,ig)
                        ze3d=xlocal8r(3,i,ig)
                        call shape8h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     elseif(lakonl(4:4).eq.'8') then
                        xi3d=xlocal8(1,i,ig)
                        et3d=xlocal8(2,i,ig)
                        ze3d=xlocal8(3,i,ig)
                        call shape8h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     elseif(lakonl(4:6).eq.'20R') then
                        xi3d=xlocal8(1,i,ig)
                        et3d=xlocal8(2,i,ig)
                        ze3d=xlocal8(3,i,ig)
                        call shape20h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     elseif(lakonl(4:4).eq.'2') then
                        xi3d=xlocal20(1,i,ig)
                        et3d=xlocal20(2,i,ig)
                        ze3d=xlocal20(3,i,ig)
                        call shape20h(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     elseif(lakonl(4:5).eq.'10') then
                        xi3d=xlocal10(1,i,ig)
                        et3d=xlocal10(2,i,ig)
                        ze3d=xlocal10(3,i,ig)
                        call shape10tet(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     elseif(lakonl(4:4).eq.'4') then
                        xi3d=xlocal4(1,i,ig)
                        et3d=xlocal4(2,i,ig)
                        ze3d=xlocal4(3,i,ig)
                        call shape4tet(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     elseif(lakonl(4:5).eq.'15') then
                        xi3d=xlocal15(1,i,ig)
                        et3d=xlocal15(2,i,ig)
                        ze3d=xlocal15(3,i,ig)
                        call shape15w(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     elseif(lakonl(4:4).eq.'6') then
                        xi3d=xlocal6(1,i,ig)
                        et3d=xlocal6(2,i,ig)
                        ze3d=xlocal6(3,i,ig)
                        call shape6w(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
                     endif
                  endif
!     
!     calculating of
!     the temperature temp
!     the static pressure pres
!     the velocity gradient vkl
!     in the integration point
!     
                  if(prlab(ii)(1:4).eq.'FLUX') then
                     temp=0.d0
                     do j1=1,3
                        vkl(0,j1)=0.d0
                     enddo
                     do i1=1,nope
                        temp=temp+shp(4,i1)*vl(0,i1)
                        do k1=1,3
                           vkl(0,k1)=vkl(0,k1)+shp(k1,i1)*vl(0,i1)
                        enddo
                     enddo
!     
!     material data (conductivity)
!     
                     call materialdata_cond(imat,ntmat_,temp,cocon,
     &                    ncocon,cond)
!     
!     determining the stress 
!     
                     dd=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
!
                     tf(0)=-cond*(vkl(0,1)*xsj2(1)+
     &                            vkl(0,2)*xsj2(2)+
     &                            vkl(0,3)*xsj2(3))
                     f(0)=f(0)+tf(0)*weight
                     tf(0)=tf(0)/dd
                     write(5,'(i10,1x,i3,1x,i3,1p,4(1x,e13.6))')nelem,
     &                    ig,i,tf(0),(coords(i1),i1=1,3)
!     
                  else
!
!                    forces and moments in sections
!                    These are obtained by integrating the stress tensor;
!                    the stresses at the integration points are interpolated
!                    from the values at the nodes of the face; these values
!                    were extrapolated from the integration point values within
!                    the element.
!
!                    interpolation of the stress tensor at the integration 
!                    point
!
                     do i1=1,3
                        do j1=i1,3
                           t(i1,j1)=0.d0
                        enddo
                     enddo
!
                     if((nope.eq.20).or.(nope.eq.8)) then
                        do i1=1,nopes
                           t(1,1)=t(1,1)+
     &                            shp2(4,i1)*stn(1,konl(ifaceq(i1,ig)))
                           t(2,2)=t(2,2)+
     &                            shp2(4,i1)*stn(2,konl(ifaceq(i1,ig)))
                           t(3,3)=t(3,3)+
     &                            shp2(4,i1)*stn(3,konl(ifaceq(i1,ig)))
                           t(1,2)=t(1,2)+
     &                            shp2(4,i1)*stn(4,konl(ifaceq(i1,ig)))
                           t(1,3)=t(1,3)+
     &                            shp2(4,i1)*stn(5,konl(ifaceq(i1,ig)))
                           t(2,3)=t(2,3)+
     &                            shp2(4,i1)*stn(6,konl(ifaceq(i1,ig)))
                         enddo
c                         write(*,*) 'printoutface'
c                         write(*,*) jface,i,t(1,1),t(2,2)
c                         write(*,*) jface,i,t(3,3),t(1,2)
c                         write(*,*) jface,i,t(1,3),t(2,3)
                     elseif((nope.eq.10).or.(nope.eq.4)) then
                        do i1=1,nopes
                           t(1,1)=t(1,1)+
     &                            shp2(4,i1)*stn(1,konl(ifacet(i1,ig)))
                           t(2,2)=t(2,2)+
     &                            shp2(4,i1)*stn(2,konl(ifacet(i1,ig)))
                           t(3,3)=t(3,3)+
     &                            shp2(4,i1)*stn(3,konl(ifacet(i1,ig)))
                           t(1,2)=t(1,2)+
     &                            shp2(4,i1)*stn(4,konl(ifacet(i1,ig)))
                           t(1,3)=t(1,3)+
     &                            shp2(4,i1)*stn(5,konl(ifacet(i1,ig)))
                           t(2,3)=t(2,3)+
     &                            shp2(4,i1)*stn(6,konl(ifacet(i1,ig)))
                        enddo
                     else
                        do i1=1,nopes
                           t(1,1)=t(1,1)+
     &                            shp2(4,i1)*stn(1,konl(ifacew(i1,ig)))
                           t(2,2)=t(2,2)+
     &                            shp2(4,i1)*stn(2,konl(ifacew(i1,ig)))
                           t(3,3)=t(3,3)+
     &                            shp2(4,i1)*stn(3,konl(ifacew(i1,ig)))
                           t(1,2)=t(1,2)+
     &                            shp2(4,i1)*stn(4,konl(ifacew(i1,ig)))
                           t(1,3)=t(1,3)+
     &                            shp2(4,i1)*stn(5,konl(ifacew(i1,ig)))
                           t(2,3)=t(2,3)+
     &                            shp2(4,i1)*stn(6,konl(ifacew(i1,ig)))
                        enddo
                     endif
                     t(2,1)=t(1,2)
                     t(3,1)=t(1,3)
                     t(3,2)=t(2,3)
!
!                    calculating the force contribution
!
                     do i1=1,3
                        df(i1)=(t(i1,1)*xsj2(1)+
     &                          t(i1,2)*xsj2(2)+
     &                          t(i1,3)*xsj2(3))*weight
                     enddo
!
!                    update total force and total moment
!
                     do i1=1,3
                        f(i1)=f(i1)+df(i1)
                     enddo
                     xm(1)=xm(1)+coords(2)*df(3)-coords(3)*df(2)
                     xm(2)=xm(2)+coords(3)*df(1)-coords(1)*df(3)
                     xm(3)=xm(3)+coords(1)*df(2)-coords(2)*df(1)
!
!                    update area, center of gravity and mean normal
!
                     dd=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
                     area=area+dd*weight
                     do i1=1,3
                        cg(i1)=cg(i1)+coords(i1)*dd*weight
                        xn(i1)=xn(i1)+xsj2(i1)*weight
                     enddo
                   endif
c                   write(*,*) 'printoutface'
c                   write(*,*) jface,i,(f(i1),i1=1,3)
c                   write(*,*) jface,i,(xm(i1),i1=1,3)
c                   write(*,*) jface,i,(cg(i1),i1=1,3)
c                   write(*,*) jface,i,(xn(i1),i1=1,3)
               enddo ! integration points in face
            enddo ! faces in surface
!
            if(prlab(ii)(1:4).eq.'FLUX') then
               write(5,*)
               write(5,123) faset(1:ipos-2),ttime+time
 123           format(' total surface flux (q) for set ',A,
     &              ' and time ',e14.7)
               write(5,*)
               write(5,'(1p,1x,e13.6)') f(0)
            else
!
!              surface set statistics
!
               write(5,*)
               write(5,130) faset(1:ipos-2),ttime+time
 130           format(' statistics for surface set ',A,
     &              ' and time ',e14.7)
!
!              total force and moment about the origin
!
               write(5,*)
               write(5,126)
!
 126           format('   total surface force (fx,fy,fz) ',
     &              'and moment about the origin(mx,my,mz)')
               write(5,*)
               write(5,'(2x,1p,6(1x,e13.6))') (f(j),j=1,3),(xm(j),j=1,3)
!
!              center of gravity and mean normal
!
               do i1=1,3
                  cg(i1)=cg(i1)/area
                  xn(i1)=xn(i1)/area
               enddo
               write(5,*)
               write(5,127)
 127           format('   center of gravity and mean normal')
               write(5,*)
               write(5,'(2x,1p,6(1x,e13.6))')(cg(j),j=1,3),(xn(j),j=1,3)
!
!              moment about the center of gravity
!
               write(5,*)
               write(5,129)
 129           format(
     &        '   moment about the center of gravity(mx,my,mz)')
               write(5,*)
               xmcg(1)=xm(1)-cg(2)*f(3)+cg(3)*f(2)
               xmcg(2)=xm(2)-cg(3)*f(1)+cg(1)*f(3)
               xmcg(3)=xm(3)-cg(1)*f(2)+cg(2)*f(1)
               write(5,'(2x,1p,6(1x,e13.6))') 
     &               xmcg(1),xmcg(2),xmcg(3)
!
!              area, normal, shear force, torque and bending moment
!
               xnormforc=f(1)*xn(1)+f(2)*xn(2)+f(3)*xn(3)
               shearforc=sqrt((f(1)-xnormforc*xn(1))**2+
     &                        (f(2)-xnormforc*xn(2))**2+
     &                        (f(3)-xnormforc*xn(3))**2)
               xtorque=xmcg(1)*xn(1)+xmcg(2)*xn(2)+xmcg(3)*xn(3)
               bendmom=sqrt((xmcg(1)-xtorque*xn(1))**2+
     &                        (xmcg(2)-xtorque*xn(2))**2+
     &                        (xmcg(3)-xtorque*xn(3))**2)
               write(5,*)
               write(5,128)
 128           format('   area, normal force (+ = tension), shear force
     &(size), torque and bending moment (size)')
               write(5,*)
               write(5,'(2x,1p,5(1x,e13.6))') area,xnormforc,shearforc,
     &              xtorque,bendmom
            endif
!     
         endif
      enddo
!     
      return
      end
