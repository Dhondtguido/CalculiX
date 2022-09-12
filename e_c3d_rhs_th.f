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
      subroutine e_c3d_rhs_th(co,nk,konl,lakonl,
     &     ff,nelem,nmethod,t0,t1,vold,nelemload,
     &     sideload,xload,nload,idist,dtime,
     &     ttime,time,istep,iinc,xloadold,reltime,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi,
     &     ielprop,prop,sti,xstateini,xstate,nstate_)
!     
!     computation of the rhs for the element with
!     the topology in konl
!     
!     ff: rhs without temperature and eigenstress contribution
!     
      implicit none
!     
      logical ivolumeforce
!
      character*8 lakonl
      character*20 sideload(*)
!
      integer konl(20),ifaceq(8,6),nelemload(2,*),nk,nelem,nmethod,
     &  nload,idist,i,j,k,i1,iflag,ipompc(*),nodempc(3,*),nmpc,
     &  jj,id,ipointer,ig,kk,nope,nopes,mint2d,ikmpc(*),ilmpc(*),
     &  mint3d,ifacet(6,4),nopev,ifacew(8,5),iinc,istep,jltyp,
     &  iscale,mi(*),ielprop(*),null,nstate_
!
      real*8 co(3,*),xl(3,20),shp(4,20),xs2(3,7),xloadold(2,*),
     &  ff(60),shpj(4,20),dxsj2,temp,press,t0(*),t1(*),coords(3),
     &  xl2(3,8),xsj2(3),shp2(7,8),vold(0:mi(2),*),xload(2,*),
     &  xi,et,ze,xsj,xsjj,t1l,ttime,time,weight,pgauss(3),tvar(2),
     &  reltime,areaj,coefmpc(*),tl2(8),prop(*),sti(6,mi(1),*),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*)
!
      real*8 dtime
!
      include "gauss.f"
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
      data null /0/
!
      tvar(1)=time
      tvar(2)=ttime+time
!
      if(lakonl(4:4).eq.'2') then
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
      else
         nope=6
         nopev=6
      endif
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
     &           null,xi,et,ze,weight)
         endif
      elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
         mint2d=4
         mint3d=8
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
      else
         mint3d=2
      endif
!
!     computation of the coordinates of the local nodes
!
      do i=1,nope
        do j=1,3
          xl(j,i)=co(j,konl(i))
        enddo
      enddo
!
!       initialisation for distributed forces
!
c      if(idist.ne.0) then
      do i=1,nope
         ff(i)=0.d0
      enddo
c     endif
!     
!     computation of the body forces
!     
      ivolumeforce=.false.
c     if(nload.gt.0) then
      call nident2(nelemload,nelem,nload,id)
      do
         if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
         if(sideload(id)(1:2).ne.'BF') then
            id=id-1
            cycle
         else
            ivolumeforce=.true.
            exit
         endif
      enddo
c     endif
!     
!     computation of the matrix: loop over the Gauss points
!     
      if(ivolumeforce) then
         do kk=1,mint3d
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
!     
!     calculation of the shape functions and their derivatives
!     in the gauss point
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
!     check the jacobian determinant
!     
            if(xsj.lt.1.d-20) then
               write(*,*) '*ERROR in e_c3d_rhs_th: nonpositive jacobian'
               write(*,*) '         determinant in element',nelem
               write(*,*)
               xsj=dabs(xsj)
               nmethod=0
            endif
!     
!     calculating the temperature in the integration
!     point
!     
            t1l=0.d0
!     
            do i1=1,nope
               t1l=t1l+shp(4,i1)*vold(0,konl(i1))
            enddo
!     
!     incorporating the jacobian determinant in the shape
!     functions
!     
            xsjj=dsqrt(xsj)
            do i1=1,nope
               shpj(1,i1)=shp(1,i1)*xsjj
               shpj(2,i1)=shp(2,i1)*xsjj
               shpj(3,i1)=shp(3,i1)*xsjj
               shpj(4,i1)=shp(4,i1)*xsj
            enddo
!     
!     computation of the right hand side
!     
!     distributed heat flux
!     
c     if(nload.gt.0) then
            call nident2(nelemload,nelem,nload,id)
            areaj=xsj*weight
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
     &                       shp(4,i1)*co(j,konl(i1))
                     enddo
                  enddo
                  jltyp=1
                  iscale=1
                  call dflux(xload(1,id),t1l,istep,iinc,tvar,
     &                 nelem,kk,pgauss,jltyp,temp,press,sideload(id),
     &                 areaj,vold,co,lakonl,konl,ipompc,nodempc,coefmpc,
     &                 nmpc,ikmpc,ilmpc,iscale,mi)
                  if((nmethod.eq.1).and.(iscale.ne.0))
     &                 xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltime
               endif
               do jj=1,nope
                  ff(jj)=ff(jj)+xload(1,id)*shpj(4,jj)*weight
               enddo
               exit
            enddo
c     endif
!     
         enddo
      endif
!     
!     distributed loads
!     
c     if(nload.eq.0) then
c     return
c     endif
      call nident2(nelemload,nelem,nload,id)
      do
         if((id.eq.0).or.(nelemload(1,id).ne.nelem)) exit
         if((sideload(id)(1:1).ne.'F').and.
     &        (sideload(id)(1:1).ne.'R').and.
     &        (sideload(id)(1:1).ne.'S')) then
            id=id-1
            cycle
         endif
         read(sideload(id)(2:2),'(i1)') ig
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
         if((nope.eq.20).or.(nope.eq.8)) then
            do i=1,nopes
               tl2(i)=vold(0,konl(ifaceq(i,ig)))
               do j=1,3
                  xl2(j,i)=co(j,konl(ifaceq(i,ig)))+
     &                 vold(j,konl(ifaceq(i,ig)))
               enddo
            enddo
         elseif((nope.eq.10).or.(nope.eq.4)) then
            do i=1,nopes
               tl2(i)=vold(0,konl(ifacet(i,ig)))
               do j=1,3
                  xl2(j,i)=co(j,konl(ifacet(i,ig)))+
     &                 vold(j,konl(ifacet(i,ig)))
               enddo
            enddo
         else
            do i=1,nopes
               tl2(i)=vold(0,konl(ifacew(i,ig)))
               do j=1,3
                  xl2(j,i)=co(j,konl(ifacew(i,ig)))+
     &                 vold(j,konl(ifacew(i,ig)))
               enddo
            enddo
         endif
!     
         do i=1,mint2d
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
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
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
               xi=gauss2d5(1,i)
               et=gauss2d5(2,i)
               weight=weight2d5(i)
            elseif((lakonl(4:4).eq.'4').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
               xi=gauss2d4(1,i)
               et=gauss2d4(2,i)
               weight=weight2d4(i)
            endif
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
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))
            areaj=dxsj2*weight
!     
            temp=0.d0
            do j=1,nopes
               temp=temp+tl2(j)*shp2(4,j)
            enddo
!     
!     for nonuniform load: determine the coordinates of the
!     point (transferred into the user subroutine)
!     
            if(sideload(id)(3:4).eq.'NU') then
               do k=1,3
                  coords(k)=0.d0
                  do j=1,nopes
                     coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
               enddo
               read(sideload(id)(2:2),'(i1)') jltyp
               jltyp=jltyp+10
               if(sideload(id)(1:1).eq.'S') then
                  iscale=1
                  call dflux(xload(1,id),temp,istep,iinc,tvar,
     &                 nelem,i,coords,jltyp,temp,press,sideload(id),
     &                 areaj,vold,co,lakonl,konl,ipompc,nodempc,
     &                 coefmpc,nmpc,ikmpc,ilmpc,iscale,mi)
                  if((nmethod.eq.1).and.(iscale.ne.0))
     &                 xload(1,id)=xloadold(1,id)+
     &                 (xload(1,id)-xloadold(1,id))*reltime
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
!     flux INTO the face is positive (input deck convention)
!     this is different from the convention in the theory
!     
                  ff(ipointer)=ff(ipointer)+shp2(4,k)*xload(1,id)
     &                 *dxsj2*weight
               elseif(sideload(id)(1:1).eq.'F') then
                  write(*,*) '*ERROR in e_c3d_rhs_th.f: no'
                  write(*,*) '       film conditions allowed'
                  write(*,*) '       in an modal dynamic calculation'
                  call exit(201)
               elseif(sideload(id)(1:1).eq.'R') then
                  write(*,*) '*ERROR in e_c3d_rhs_th.f: no'
                  write(*,*) '       radiation conditions allowed'
                  write(*,*) '       in an modal dynamic calculation'
                  call exit(201)
               endif
            enddo
!     
         enddo
!     
         id=id-1
      enddo
!     
      return
      end
      


