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
      subroutine frictionheating(ne0,ne,ipkon,lakon,ielmat,mi,elcon,
     &  ncmat_,ntmat_,kon,islavsurf,pmastsurf,springarea,co,vold,
     &  veold,pslavsurf,xloadact,nload,nload_,nelemload,iamload,
     &  idefload,sideload,stx,nam,time,ttime,matname,istep,iinc)
!
!     determines the effect of friction heating
!
      implicit none
!
      character*8 lakon(*),lakonl
      character*20 label,sideload(*)
      character*80 matname(*)
!
      integer i,j,k,ne0,ne,indexe,ipkon(*),imat,mi(*),ielmat(mi(3),*),
     &  ncmat_,ntmat_,kon(*),nope,igauss,jfaces,ifaces,nelems,ifacem,
     &  nelemm,jfacem,islavsurf(2,*),nopes,nopem,konl(20),iflag,
     &  mint2d,iamplitude,isector,nload,nload_,nelemload(2,*),nam,
     &  iamload(2,*),idefload(*),istep,iinc
!
      real*8 elcon(0:ncmat_,ntmat_,*),pressure,stx(6,mi(1),*),
     &  pmastsurf(6,*),area,springarea(2,*),pl(3,20),co(3,*),
     &  vold(0:mi(2),*),areaslav,xi,et,vels(3),veold(0:mi(2),*),
     &  xsj2m(3),xs2m(3,7),shp2m(7,9),xsj2s(3),xs2s(3,7),shp2s(7,9),
     &  areamast,pslavsurf(3,*),value,velm(3),um,xloadact(2,*),weight,
     &  shear,vnorm,f,eta,timeend(2),time,ttime,coords(3),xl(3,20)
!
!
!
      include "gauss.f"
!
!     mortar=1 is assumed (face-to-face penalty contact)
!     ithermal=3 is assumed
!
      iamplitude=0
      isector=0
!
      do i=ne0+1,ne
         imat=ielmat(1,i)
!
!        heat conversion factor
!
         eta=elcon(9,1,imat)
!
!        surface weighting factor
!
         f=elcon(10,1,imat)
!
!        velocity
!
         vnorm=elcon(11,1,imat)
!
!        friction coefficient
!
         um=elcon(6,1,imat)
!
         pressure=stx(4,1,i)
         if(pressure.lt.0.d0) cycle
!
         shear=dsqrt(stx(5,1,i)**2+stx(6,1,i)**2)
         if(vnorm.lt.-0.5d0) then
!
!           if ||v||<0 => take differential velocity from the results
!           no heat generation if no slip
!
            if(shear.lt.um*pressure*0.95d0) cycle
         endif
!
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         nope=kon(ipkon(i))
         nopem=ichar(lakonl(8:8))-48
         nopes=nope-nopem
!
         igauss=kon(indexe+nope+1)
         jfaces=kon(indexe+nope+2)
!
!        slave face
!
         ifaces=islavsurf(1,jfaces)
         nelems=int(ifaces/10.d0)
         jfaces=ifaces-10*nelems
!
!        master face
!
         ifacem=int(pmastsurf(3,igauss))
         nelemm=int(ifacem/10.d0)
         jfacem=ifacem-10*nelemm
!
!        contact area
!
         area=springarea(1,igauss)
!
!        slave and master nodes
!
         do j=1,nope
            konl(j)=kon(indexe+j)
            do k=1,3
               pl(k,j)=co(k,konl(j))+vold(k,konl(j))
            enddo
         enddo
!
!        user subroutine called if vnorm=-0.01d0
!
         if((vnorm.lt.0.d0).and.(vnorm.gt.-0.5d0)) then
!
            xi=pslavsurf(1,igauss)
            et=pslavsurf(2,igauss)
!
            do j=nopem+1,nopem+nopes
               konl(j)=kon(indexe+j)
               do k=1,3
                  xl(k,j)=co(k,konl(j))
               enddo
            enddo
!
!           determining the jacobian vector on the surface 
!
            iflag=1
            if(nopes.eq.8) then
               call shape8q(xi,et,xl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
            elseif(nopes.eq.4) then
               call shape4q(xi,et,xl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
            elseif(nopes.eq.6) then
               call shape6tri(xi,et,xl(1,nopem+1),xsj2s,xs2s,shp2s,
     &                        iflag)
            else
               call shape3tri(xi,et,xl(1,nopem+1),xsj2s,xs2s,shp2s,
     &                        iflag)
            endif
!
!           position of the slave integration point
!
            do j=1,3
               coords(j)=0.d0
               do k=1,nopes
                  coords(j)=coords(j)+shp2s(4,k)*xl(j,nopem+k)
               enddo
            enddo
!
            timeend(1)=time
            timeend(2)=ttime+time
            call fricheat(eta,f,vnorm,timeend,matname(imat),i,
     &                    nelems,jfaces,nelemm,jfacem,um,
     &                    istep,iinc,area,pressure,coords)
         endif
!
!        heat flux into the slave face
!
         if(nopes.eq.8) then
            mint2d=9
         elseif(nopes.eq.6) then
            mint2d=3
         elseif(nopes.eq.4) then
            mint2d=4
         else
            mint2d=1
         endif
!
!        calculating the area of the slave face
!
         areaslav=0.d0
!
         do j=1,mint2d
            if(nopes.eq.8) then
               xi=gauss2d3(1,j)
               et=gauss2d3(2,j)
               weight=weight2d3(j)
            elseif(nopes.eq.6) then
               xi=gauss2d5(1,j)
               et=gauss2d5(2,j)
               weight=weight2d5(j)
            elseif(nopes.eq.4) then
               xi=gauss2d2(1,j)
               et=gauss2d2(2,j)
               weight=weight2d2(j)
            else
               xi=gauss2d4(1,j)
               et=gauss2d4(2,j)
               weight=weight2d4(j)
            endif
!
            iflag=2
            if(nopes.eq.8) then
               call shape8q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
            elseif(nopes.eq.4) then
               call shape4q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
            elseif(nopes.eq.6) then
               call shape6tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &iflag)
            else
               call shape3tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &iflag)
            endif
!
            areaslav=areaslav+dsqrt(xsj2s(1)**2+
     &                              xsj2s(2)**2+
     &                              xsj2s(3)**2)
         enddo
!
         label(1:20)='S                   '
         write(label(2:2),'(i1)') jfaces
!
         if(vnorm.gt.0.d0) then
            value=um*pressure*vnorm*eta*f*area/areaslav
            call loadadd(nelems,label,value,nelemload,sideload,xloadact,
     &                nload,nload_,iamload,iamplitude,nam,isector,
     &                idefload)
         elseif(vnorm.lt.-0.5d0) then
!
!           calculate the differential velocity
!
!           determining the slave velocity
!
            xi=pslavsurf(1,igauss)
            et=pslavsurf(2,igauss)
            iflag=1
            if(nopes.eq.8) then
               call shape8q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
            elseif(nopes.eq.4) then
               call shape4q(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,iflag)
            elseif(nopes.eq.6) then
               call shape6tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &                        iflag)
            else
               call shape3tri(xi,et,pl(1,nopem+1),xsj2s,xs2s,shp2s,
     &                        iflag)
            endif
!
            do k=1,3
               vels(k)=0.d0
               do j=1,nopes
                  vels(k)=vels(k)+shp2s(4,j)*
     &                            veold(k,konl(nopem+j))
               enddo
            enddo
!
!           determining the master velocity
!
            xi=pmastsurf(1,igauss)
            et=pmastsurf(2,igauss)
            iflag=1
            if(nopem.eq.8) then
               call shape8q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            elseif(nopem.eq.4) then
               call shape4q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            elseif(nopem.eq.6) then
               call shape6tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            else
               call shape3tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            endif
!
            do k=1,3
               velm(k)=0.d0
               do j=1,nopem
                  velm(k)=velm(k)+shp2m(4,j)*
     &                            veold(k,konl(j))
               enddo
            enddo
!
            vnorm=dsqrt((vels(1)-velm(1))**2+
     &                  (vels(2)-velm(2))**2+
     &                  (vels(3)-velm(3))**2)
            value=um*pressure*vnorm*eta*f*area/areaslav
            call loadadd(nelems,label,value,nelemload,sideload,xloadact,
     &                nload,nload_,iamload,iamplitude,nam,isector,
     &                idefload)
         endif
!
!        heat flux into the master face
!
         if(nopem.eq.8) then
            mint2d=9
         elseif(nopem.eq.6) then
            mint2d=3
         elseif(nopem.eq.4) then
            mint2d=4
         else
            mint2d=1
         endif
!
!        calculating the area of the slave face
!
         areamast=0.d0
!
         do j=1,mint2d
            if(nopem.eq.8) then
               xi=gauss2d3(1,j)
               et=gauss2d3(2,j)
               weight=weight2d3(j)
            elseif(nopem.eq.6) then
               xi=gauss2d5(1,j)
               et=gauss2d5(2,j)
               weight=weight2d5(j)
            elseif(nopem.eq.4) then
               xi=gauss2d2(1,j)
               et=gauss2d2(2,j)
               weight=weight2d2(j)
            else
               xi=gauss2d4(1,j)
               et=gauss2d4(2,j)
               weight=weight2d4(j)
            endif
!
            iflag=2
            if(nopem.eq.8) then
               call shape8q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            elseif(nopem.eq.4) then
               call shape4q(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            elseif(nopem.eq.6) then
               call shape6tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            else
               call shape3tri(xi,et,pl,xsj2m,xs2m,shp2m,iflag)
            endif
!
            areamast=areamast+dsqrt(xsj2m(1)**2+
     &                              xsj2m(2)**2+
     &                              xsj2m(3)**2)
         enddo
!
         label(1:20)='S                   '
         write(label(2:2),'(i1)') jfacem
!
!        at this point vnorm was either given by the user or
!        calculated for the slave surface (differential velocity)
!
         value=um*pressure*vnorm*eta*(1.d0-f)*area/areamast
         call loadadd(nelemm,label,value,nelemload,sideload,xloadact,
     &                nload,nload_,iamload,iamplitude,nam,isector,
     &                idefload)
      enddo
!
      return
      end
      
