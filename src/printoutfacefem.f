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
      subroutine printoutfacefem(co,ntmat_,vold,shcon,
     &     nshcon,compressible,istartset,iendset,ipkon,
     &     lakon,kon,ialset,prset,ttime,nset,set,nprint,prlab,ielmat,mi,
     &     nelnew)
!     
!     calculation and printout of the lift and drag forces
!     
      implicit none
!     
      integer compressible
!     
      character*8 lakonl,lakon(*)
      character*6 prlab(*)
      character*80 faset
      character*81 set(*),prset(*)
!     
      integer konl(8),ifaceq(8,6),nelem,ii,nprint,i,j,i1,i2,j1,id,
     &     k1,jj,ig,nshcon(*),ntmat_,nope,nopes,nelef,nelnew(*),
     &     imat,mint2d,ifacet(6,4),ifacew(8,5),iflag,indexe,jface,
     &     istartset(*),iendset(*),ipkon(*),kon(*),iset,ialset(*),nset,
     &     ipos,mi(*),ielmat(mi(3),*)
!     
      real*8 co(3,*),xl(3,8),shp(4,8),xs2(3,7),dvi,f(3),
     &     vkl(3,3),t(3,3),div,shcon(0:3,ntmat_,*),
     &     voldl(0:mi(2),8),xl2(3,8),xsj2(3),
     &     shp2(7,8),vold(0:mi(2),*),xi,et,xsj,temp,xi3d,et3d,ze3d,
     &     weight,xlocal20(3,9,6),xlocal4(3,1,4),xlocal10(3,3,4),
     &     xlocal6(3,1,5),xlocal15(3,4,5),xlocal8(3,4,6),
     &     xlocal8r(3,1,6),ttime,pres,tf(3),tn,tt,dd,coords(3)
!     
      include "gauss.f"
      include "xlocal.f"
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
      data iflag /3/
!     
!     initialisierung forces
!     
      do i=1,3
        f(i)=0.d0
      enddo
!     
      do ii=1,nprint
!     
!     total drag
!     
        if(prlab(ii)(1:4).eq.'DRAG') then
!     
          ipos=index(prset(ii),' ')
          faset='                    '
          faset(1:ipos-1)=prset(ii)(1:ipos-1)
!     
!     printing the header
!     
          write(5,*)
          write(5,120) faset(1:ipos-2),ttime
 120      format(' surface stress vector (tx,ty,tz), normal stress, shea
     &r stress and coordinates for set ',A,' and time ',e14.7)
          write(5,*)
!     
!     printing the data
!     
c          do iset=1,nset
c            if(set(iset).eq.prset(ii)) exit
c          enddo
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
            nelef=nelnew(nelem)
!            
            lakonl=lakon(nelef)
            indexe=ipkon(nelef)
            imat=ielmat(1,nelef)
!     
            if(lakonl(4:4).eq.'8') then
              nope=8
              nopes=4
            elseif(lakonl(4:4).eq.'4') then
              nope=4
              nopes=3
            elseif(lakonl(4:4).eq.'6') then
              nope=6
            endif
!     
            if(lakonl(4:5).eq.'8R') then
              mint2d=1
            elseif(lakonl(4:4).eq.'8') then
              mint2d=4
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
!     temperature, velocity and auxiliary variables
!     (rho*energy density, rho*velocity and rho)
!     
            do i1=1,nope
              do i2=0,4
                voldl(i2,i1)=vold(i2,konl(i1))
              enddo
            enddo
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
!     
            if(nope.eq.8) then
              do i=1,nopes
                do j=1,3
                  xl2(j,i)=co(j,konl(ifaceq(i,ig)))
                enddo
              enddo
            elseif(nope.eq.4) then
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
!     
            do i=1,mint2d
!     
!     local coordinates of the surface integration
!     point within the surface local coordinate system
!     
              if((lakonl(4:5).eq.'8R').or.
     &             ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
                xi=gauss2d1(1,i)
                et=gauss2d1(2,i)
                weight=weight2d1(i)
              elseif(lakonl(4:4).eq.'8') then
                xi=gauss2d2(1,i)
                et=gauss2d2(2,i)
                weight=weight2d2(i)
              elseif((lakonl(4:4).eq.'4').or.
     &               ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
                xi=gauss2d4(1,i)
                et=gauss2d4(2,i)
                weight=weight2d4(i)
              endif
!     
!     local surface normal
!     
              if(nopes.eq.4) then
                call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
              else
                call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
              endif
!     
!     global coordinates of the integration point
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
              elseif(lakonl(4:4).eq.'4') then
                xi3d=xlocal4(1,i,ig)
                et3d=xlocal4(2,i,ig)
                ze3d=xlocal4(3,i,ig)
                call shape4tet(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
              elseif(lakonl(4:4).eq.'6') then
                xi3d=xlocal6(1,i,ig)
                et3d=xlocal6(2,i,ig)
                ze3d=xlocal6(3,i,ig)
                call shape6w(xi3d,et3d,ze3d,xl,xsj,shp,iflag)
              endif
!     
!     calculating of
!     the temperature temp
!     the static pressure pres
!     the velocity gradient vkl
!     in the integration point
!     
              temp=0.d0
              pres=0.d0
              do i1=1,3
                do j1=1,3
                  vkl(i1,j1)=0.d0
                enddo
              enddo
              do i1=1,nope
                temp=temp+shp(4,i1)*voldl(0,i1)
                pres=pres+shp(4,i1)*voldl(4,i1)
                do j1=1,3
                  do k1=1,3
                    vkl(j1,k1)=vkl(j1,k1)+shp(k1,i1)*voldl(j1,i1)
                  enddo
                enddo
              enddo
              if(compressible.eq.1) div=vkl(1,1)+vkl(2,2)+vkl(3,3)
!     
!     material data (density, dynamic viscosity, heat capacity and
!     conductivity)
!     
              call materialdata_dvifem(imat,ntmat_,temp,shcon,nshcon,
     &             dvi)
!     
!     determining the stress 
!     
              do i1=1,3
                do j1=1,3
                  t(i1,j1)=vkl(i1,j1)+vkl(j1,i1)
                enddo
                if(compressible.eq.1) 
     &               t(i1,i1)=t(i1,i1)-2.d0*div/3.d0
              enddo
!     
              dd=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &             xsj2(3)*xsj2(3))
              do i1=1,3
                tf(i1)=dvi*(t(i1,1)*xsj2(1)+t(i1,2)*xsj2(2)+
     &               t(i1,3)*xsj2(3))-pres*xsj2(i1)
                f(i1)=f(i1)+tf(i1)*weight
                tf(i1)=tf(i1)/dd
              enddo
              tn=(tf(1)*xsj2(1)+tf(2)*xsj2(2)+tf(3)*xsj2(3))/dd
              tt=dsqrt((tf(1)-tn*xsj2(1)/dd)**2+
     &             (tf(2)-tn*xsj2(2)/dd)**2+
     &             (tf(3)-tn*xsj2(3)/dd)**2)
              write(5,'(i6,1x,i3,1x,i3,1p,8(1x,e11.4))') nelem,ig,i,
     &             (tf(i1),i1=1,3),tn,tt,(coords(i1),i1=1,3)
!     
            enddo
          enddo
!     
          write(5,*)
          write(5,121) faset(1:ipos-2),ttime
 121      format(' total surface force (fx,fy,fz) for set ',A,
     &         ' and time ',e14.7)
          write(5,*)
          write(5,'(1p,3(1x,e11.4))') (f(j),j=1,3)
!     
        endif
      enddo
!     
      return
      end
      
      
