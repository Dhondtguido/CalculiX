!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine advecforc(nope,voldl,ithermal,xl,nelemload,nelemadvec,
     &  nload,lakon,xload,istep,time,ttime,dtime,sideload,vold,mi,
     &  xloadold,reltime,nmethod,tnl,iinc,iponoel,inoel,ielprop,prop,
     &  ielmat,shcon,nshcon,rhcon,nrhcon,ntmat_,ipkon,kon,cocon,ncocon,
     &  ipobody,ibody,xbody)
!
!     calculates the stiffness of an advective element.
!     An advective element consists of a face with a forced convection
!     film condition and a network node
!
      implicit none
!
      character*8 lakonl,lakon(*) 
      character*20 sideload(*),sideloadl
!
      integer nope,i,ithermal(*),j,nelemload(2,*),nelemadvec,nload,id,
     &  nelem,ig,mint2d,iflag,istep,jltyp,nfield,mi(*),nmethod,k,iinc,
     &  node,nopes,iponoel(*),inoel(2,*),ielprop(*),ielmat(mi(3),*),
     &  ipkon(*),
     &  nshcon(*),nrhcon(*),ntmat_,kon(*),ncocon(2,*),ipobody(2,*),
     &  ibody(3,*)
!
      real*8 tl2(9),voldl(0:mi(2),20),xl(3,9),sinktemp,xi,et,weight,
     &  xl2(3,8),xsj2(3),shp2(7,9),coords(3),xs2(3,7),dxsj2,areaj,
     &  temp,xload(2,*),timeend(2),time,ttime,dtime,field,reltime,
     &  vold(0:mi(2),*),xloadold(2,*),tnl(9),tnlref,prop(*),
     &  shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),cocon(0:6,ntmat_,*),
     &  xbody(7,*),heatnod,heatfac
!
      include "gauss.f"
!
      data iflag /2/
!
      timeend(1)=time
      timeend(2)=ttime+time
!
      heatfac=0.d0
!
!     number of nodes in the advective face
!
      nopes=nope-1
!
!     temperature and displacements in the element's nodes
!
      do i=1,nope
         tl2(i)=voldl(0,i)
      enddo
      do i=1,nopes
         if(ithermal(2).eq.2) then
            do j=1,3
               xl2(j,i)=xl(j,i)
            enddo
         else
            do j=1,3
               xl2(j,i)=xl(j,i)+voldl(j,i)
            enddo
         endif
      enddo
!
      call nident2(nelemload,nelemadvec,nload,id)
!
!     the second entry in nelemload points to the original
!     film loading
!
      id=nelemload(2,id)
!
!     number of the original element
!
      nelem=nelemload(1,id)
      lakonl=lakon(nelem)
!
!     the next line originally ran:
!     read(sideload(id)(2:2),'(i1)') ig
!     it was replaced since the read statement might
!     cause problems in the multithreading version
!
      ig=ichar(sideload(id)(2:2))-48
!
!     number of integration points
!
      if(lakonl(4:5).eq.'8R') then
         mint2d=1
      elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')) then
         if((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'E')) then
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
      elseif(lakonl(4:5).eq.'15') then
         if(ig.le.2) then
            mint2d=3
         else
            mint2d=4
         endif
      elseif(lakonl(4:4).eq.'6') then
         mint2d=1
      endif
!
      do i=1,mint2d
!
!     copying the sink temperature to ensure the same
!     value in each integration point (sinktemp can be
!     changed in subroutine film: requirement from the
!     thermal people)
!     
         sinktemp=tl2(nope)
!     
         if((lakonl(4:5).eq.'8R').or.
     &        ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
            xi=gauss2d1(1,i)
            et=gauss2d1(2,i)
            weight=weight2d1(i)
         elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R').or.
     &           ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
            xi=gauss2d2(1,i)
            et=gauss2d2(2,i)
            weight=weight2d2(i)
         elseif(lakonl(4:4).eq.'2') then
            xi=gauss2d3(1,i)
            et=gauss2d3(2,i)
            weight=weight2d3(i)
         elseif((lakonl(4:5).eq.'10').or.
     &           ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
            xi=gauss2d5(1,i)
            et=gauss2d5(2,i)
            weight=weight2d5(i)
         elseif((lakonl(4:4).eq.'4').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
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
     &        xsj2(3)*xsj2(3))
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
         if((sideload(id)(3:4).eq.'NU').or.
     &        (sideload(id)(5:6).eq.'NU')) then
            do k=1,3
               coords(k)=0.d0
               do j=1,nopes
                  coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
               enddo
            enddo
            jltyp=ichar(sideload(id)(2:2))-38
            node=nelemload(2,id)
            sideloadl=sideload(id)
            sideloadl(1:1)='F'
            call film(xload(1,id),sinktemp,temp,istep,
     &           iinc,timeend,nelem,i,coords,jltyp,field,nfield,
     &           sideloadl,node,areaj,vold,mi,
     &           ipkon,kon,lakon,iponoel,inoel,ielprop,prop,ielmat,
     &           shcon,nshcon,rhcon,nrhcon,ntmat_,cocon,ncocon,
     &           ipobody,xbody,ibody,heatnod,heatfac)
         endif
!
         tnlref=(xload(1,id)*(sinktemp-temp)+heatfac)*areaj
         do j=1,nopes
            tnl(j)=tnl(j)-tnlref*shp2(4,j)
         enddo
         tnl(nope)=tnl(nope)+tnlref
      enddo
!
!     for some axisymmetric and plane strain elements only half the
!     number of integration points was used
!
      if(((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R')).and.
     &    ((lakonl(7:7).eq.'A').or.(lakonl(7:7).eq.'E'))) then
         do i=1,nope
            tnl(i)=2.d0*tnl(i)
         enddo
      endif
!     
      return
      end
      
